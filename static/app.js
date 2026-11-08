// Find an interface element by its unique HTML identifier.
function byId(id) {
  return document.getElementById(id);
}

// Escape source text before inserting it into HTML or SVG markup.
function escapeHtml(value) {
  return String(value ?? "").replace(/[&<>"']/g, escapeCharacter);
}

// Replace one HTML-sensitive character with its safe entity representation.
function escapeCharacter(character) {
  const entities = {
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
    "'": "&#39;",
  };
  return entities[character];
}

let neighboringGenes = [];
let neighborOptions = [];
let activeNeighborOption = -1;
let variantSuggestionRequest = 0;
let variantSuggestionTimer;
let variantOptions = [];
let activeVariantOption = -1;
let regionQuery = "";
let expressionData = null;
let expressionExport = null;
if (
  location.pathname === "/" &&
  new URLSearchParams(location.hash.slice(1)).has("gene") &&
  new URLSearchParams(location.hash.slice(1)).get("page") !== "home"
)
  history.replaceState(null, "", "/gwas" + location.hash);
let selectedGene = "";
let currentPage = 1;
let totalAssociations = 0;
let latestRequest = 0;
let expressionRequest = 0;
let expressionCatalogLoaded = false;
let initialExpressionDataset = new URLSearchParams(location.hash.slice(1)).get(
  "expression",
);

// Fetch a JSON response and report readable API errors to the caller.
async function fetchJson(path) {
  const response = await fetch(path);
  const data = await response.json();
  if (!response.ok) {
    throw new Error(data.error || "Request failed");
  }
  return data;
}

// Display a progress message or an error in the accessible status area.
function setStatus(text, isError = false) {
  byId("message").textContent = text;
  byId("message").className = isError ? "error" : "";
}

// Search recognized gene identifiers and ask the user to resolve ambiguous matches.
async function searchGenes() {
  const requestId = ++latestRequest;
  setStatus("Searching…");
  byId("result").hidden = true;
  byId("matches").replaceChildren();

  try {
    const query = encodeURIComponent(byId("query").value);
    const data = await fetchJson(`/api/search?q=${query}`);
    if (requestId !== latestRequest) return;

    if (!data.matches.length) {
      setStatus(
        "No matching identifier in this annotation. Try a RAP or mapped MSU ID.",
      );
      return;
    }
    if (data.matches.length === 1) {
      await loadRegion(data.matches[0].id, 1);
      return;
    }

    const limitNotice = data.truncated ? " (first 50; refine your search)" : "";
    setStatus(
      `${data.matches.length} matching genes${limitNotice}. Symbols can be shared aliases; compare IDs and descriptions before selecting.`,
    );
    for (const gene of data.matches) {
      const button = document.createElement("button");
      const matched = gene.matched_identifier
        ? `Matched ${gene.matched_identifier} (${gene.matched_identifier_kind}); `
        : "";
      button.className = "gene-choice";
      button.innerHTML = `<strong>${escapeHtml(gene.id)} · ${escapeHtml(gene.symbol || "unnamed")} · chr${gene.chromosome}</strong><span>${escapeHtml(matched)}${escapeHtml(gene.description)}</span>`;
      button.dataset.gene = gene.id;
      button.onclick = chooseGene;
      byId("matches").append(button);
    }
  } catch (error) {
    if (requestId === latestRequest) setStatus(error.message, true);
  }
}

// Render one recognized gene alias together with its identifier type.
function renderAlias(alias) {
  return `${escapeHtml(alias.name)} (${escapeHtml(alias.kind)})`;
}

// Render the selected gene's coordinates, description and recognized aliases.
function renderGene(gene, aliases) {
  const symbol = gene.symbol ? `· ${escapeHtml(gene.symbol)}` : "";
  return `
    <span class="muted">SELECTED GENE · ${escapeHtml(gene.assembly_id)}</span>
    <div class="id">${escapeHtml(gene.id)} ${symbol}</div>
    <p>${escapeHtml(gene.description)}</p>
    <p class="muted">
      Chr ${gene.chromosome}: ${gene.start.toLocaleString()}–${gene.end.toLocaleString()}
      · strand ${gene.strand} · RAP-DB 2026-02-05
    </p>
    <p class="gene-tools"><a href="/variation" data-route="variation">View variants</a> · <a href="/gwas" data-route="gwas">View associations</a> · <a href="/expression" data-route="expression">View expression</a></p>
    <details>
      <summary>Recognized identifiers and aliases</summary>
      ${aliases.map(renderAlias).join("<br>")}
    </details>
  `;
}

// Render one published association without substituting values for missing statistics.
function renderAssociation(association) {
  const reviewNotice = association.curation_note
    ? `<small class="curation-warning">${escapeHtml(association.curation_note)}</small>`
    : `<small>${escapeHtml(association.review_status)}</small>`;
  const allele = association.effect_allele
    ? `Effect allele: ${escapeHtml(association.effect_allele)}`
    : "Effect allele not reported";
  const frequency = association.effect_allele_frequency
    ? `Frequency: ${escapeHtml(association.effect_allele_frequency)}`
    : "Frequency not reported";
  const repeatNotice =
    association.source_record_count > 1
      ? `${association.source_record_count} identical evidence rows grouped`
      : "One source evidence row";
  const effect =
    association.effect_size == null
      ? "Effect size not reported"
      : `Effect size: ${escapeHtml(association.effect_size)}`;
  const otherAllele =
    association.other_allele == null
      ? "Other allele not reported"
      : `Other allele: ${escapeHtml(association.other_allele)}`;
  const tableAlleles = association.original_minor_allele
    ? `<br>Original Table 1: major ${escapeHtml(association.original_major_allele)} / minor ${escapeHtml(association.original_minor_allele)}; MAF ${escapeHtml(association.original_minor_allele_frequency)} (effect orientation unverified)`
    : "";
  const variance = association.r2_percent
    ? `<br>Reported R²: ${escapeHtml(association.r2_percent)}%`
    : "";
  const publication = association.publication_url
    ? `<a href="${escapeHtml(association.publication_url)}" target="_blank" rel="noopener">
         ${escapeHtml(association.title)}
       </a>`
    : "Publication identifier unavailable";

  return `
    <tr>
      <td>${escapeHtml(association.trait)}<small>${escapeHtml(association.ontology_id)}</small>${reviewNotice}</td>
      <td>
        chr${association.chromosome}:${association.position.toLocaleString()}
        <small>Source row ${association.source_row} · record ${escapeHtml(association.id)}<br>${repeatNotice}</small>
      </td>
      <td>
        ${escapeHtml(association.relation)}
        <small>${association.distance_bp.toLocaleString()} bp from gene span</small>
      </td>
      <td>${escapeHtml(association.p_value_text)}<small>${escapeHtml(association.statistic_kind)}</small></td>
      <td>${allele}<small>${frequency}<br>${otherAllele}<br>${effect}${variance}${tableAlleles}</small></td>
      <td>
        ${escapeHtml(association.source_study_id)} · ${escapeHtml(association.model)}
        <small>
          N = ${escapeHtml(association.sample_size)} · ${escapeHtml(association.environment)}
          · ${escapeHtml(association.tissue)}
        </small>
        ${publication}
        <details>
          <summary>Coordinate provenance</summary>
          Display assembly: ${escapeHtml(association.assembly)}<br>
          ${escapeHtml(association.coordinate_status)}<br>
          Atlas reference: ${escapeHtml(association.source_reference)}<br>
          Original reported location: ${escapeHtml(association.reported_location)}<br>
          <a href="${escapeHtml(association.source_url)}" target="_blank" rel="noopener">Source dataset</a> · <a href="${escapeHtml(association.source_record_url)}" target="_blank" rel="noopener">Original source row</a>
        </details>
      </td>
    </tr>
  `;
}

// Render one neighboring gene as a selectable row with its reference orientation.
function renderNeighbor(gene) {
  const selectedLabel = gene.id === selectedGene ? " (selected)" : "";
  return `
    <tr>
      <td><button class="link" data-gene="${escapeHtml(gene.id)}">${escapeHtml(gene.id)}${selectedLabel}</button></td>
      <td>${escapeHtml(gene.symbol || "—")}</td>
      <td>${gene.start.toLocaleString()}–${gene.end.toLocaleString()} (${gene.strand})</td>
      <td>${escapeHtml(gene.description)}</td>
    </tr>
  `;
}

// Load a gene's regional evidence while ignoring responses from superseded requests.
async function loadRegion(geneId, page = 1) {
  closeVariantSuggestions();
  closeNeighborSuggestions();
  const requestId = ++latestRequest;
  setStatus("Loading gene and regional evidence…");

  try {
    if (geneId !== selectedGene) {
      byId("variant-search").value = "";
      closeVariantSuggestions();
      byId("study").value = "";
      byId("publication").value = "";
    }
    const query = new URLSearchParams({
      gene: geneId,
      window: byId("window").value,
      trait: byId("trait").value,
      page,
      study: byId("study").value,
      publication: byId("publication").value,
      maximum_p: byId("maximum-p").value,
      position: byId("position").value,
      include_repeats: byId("include-repeats").checked,
      variant_query:
        document.body.dataset.page === "variation" ? byId("variant-search").value : "",
    });
    regionQuery = query.toString();
    const data = await fetchJson(`/api/region?${query}`);
    if (requestId !== latestRequest) return;

    selectedGene = geneId;
    currentPage = data.page;
    updateStudyOptions(data);
    totalAssociations = data.total_associations;
    byId("matches").replaceChildren();
    byId("result").hidden = false;
    setStatus("");
    byId("gene").innerHTML = renderGene(data.gene, data.aliases);
    for (const link of byId("gene").querySelectorAll("[data-route]"))
      link.onclick = navigateTool;

    const plotNotice =
      totalAssociations > data.plot_count
        ? " strongest records; table and CSV include all matching records"
        : "";
    byId("region").textContent =
      `Chr ${data.gene.chromosome}: ${data.start.toLocaleString()}–${data.end.toLocaleString()} · ` +
      `${data.neighbors.length} genes · ${totalAssociations.toLocaleString()} association records. ` +
      `Plot shows ${data.plot_count.toLocaleString()}${plotNotice}.`;

    const traitNotice = data.trait_filter ? ` matching “${data.trait_filter}”` : "";
    byId("counts").textContent =
      `${totalAssociations.toLocaleString()} evidence records${traitNotice} from ${data.source_record_total.toLocaleString()} matching source rows; ${data.repeated_rows_hidden} exact repeats grouped. ` +
      "Records may share SNP positions, traits or study samples; counts are not independent replications.";
    byId("associations").innerHTML = data.associations.length
      ? data.associations.map(renderAssociation).join("")
      : '<tr><td colspan="6">No records in the imported dataset for this window and trait filter. This does not establish absence of association.</td></tr>';
    byId("variation-counts").textContent =
      `${data.total_associations.toLocaleString()} matching observations; showing ${data.associations.length} on this page. Use the page controls below to explore more.`;
    byId("variants").innerHTML = data.associations.length
      ? data.associations.map(renderVariant).join("")
      : '<tr><td colspan="5">No recorded variants match these filters. Try a wider window or reset the filters; this does not imply an absence of variation.</td></tr>';
    neighboringGenes = data.neighbors;
    byId("neighbors").innerHTML = data.neighbors.map(renderNeighbor).join("");
    bindGeneControls(byId("neighbors"));
    filterNeighborTable();

    byId("prev").disabled = currentPage === 1;
    byId("next").disabled = currentPage * data.page_size >= totalAssociations;
    byId("page").textContent =
      `Page ${currentPage} of ${Math.max(1, Math.ceil(totalAssociations / data.page_size))}`;
    byId("page-number").value = currentPage;
    byId("page-number").max = data.page_count;
    byId("export").href = `/api/export?${query}`;
    drawNeighborhood(data);
    renderExpressionControls(data.neighbors);
    loadExpression();
    history.replaceState(
      null,
      "",
      "#" +
        new URLSearchParams({
          gene: geneId,
          window: data.window,
          expression: byId("expression-dataset").value,
        }),
    );
    updateToolLinks();
  } catch (error) {
    if (requestId === latestRequest) setStatus(error.message, true);
  }
}

// Populate region-specific study and publication choices without discarding active selections.
function updateStudyOptions(data) {
  const study = data.filters.study;
  const publication = data.filters.publication;
  byId("study").replaceChildren(new Option("All studies", ""));
  for (const item of data.study_options) {
    byId("study").add(new Option(`${item.source_study_id} · ${item.model}`, item.id));
  }
  byId("publication").replaceChildren(new Option("All publications", ""));
  for (const item of data.publication_options) {
    byId("publication").add(
      new Option(`${item.pmid} · ${item.title.slice(0, 90)}`, item.pmid),
    );
  }
  byId("study").value = study;
  byId("publication").value = publication;
}

// Clear evidence filters while preserving the selected gene and genomic window.
function resetFilters() {
  byId("variant-search").value = "";
  closeVariantSuggestions();
  byId("trait").value = "";
  byId("study").value = "";
  byId("publication").value = "";
  byId("maximum-p").value = "1";
  byId("position").value = "";
  byId("include-repeats").checked = false;
  applyFilters();
}

// Navigate to the requested evidence page within the available page range.
function goToPage() {
  const page = Number(byId("page-number").value);
  if (!Number.isInteger(page) || page < 1) {
    setStatus("Enter a positive whole page number.", true);
    return;
  }
  loadRegion(selectedGene, page);
}

// Find a chart lane with enough horizontal space for the next gene label.
function findAvailableLane(laneEnds, left) {
  for (let index = 0; index < laneEnds.length; index++) {
    if (laneEnds[index] < left - 145) return index;
  }
  return laneEnds.length;
}

// Draw association strengths above selectable gene arrows using reference coordinates.
function drawNeighborhood(data) {
  const svg = byId("chart");
  const width = 1100;
  const leftMargin = 65;
  const rightMargin = 30;
  const topMargin = 25;
  const plotHeight = 160;
  let maximum = 1;
  for (const association of data.plot_associations) {
    maximum = Math.max(maximum, association.neg_log10_p);
  }

  // Convert a genomic position into a horizontal chart coordinate.
  function x(position) {
    const fraction = (position - data.start) / (data.end - data.start);
    return leftMargin + fraction * (width - leftMargin - rightMargin);
  }

  // Convert an association strength into a vertical chart coordinate.
  function y(value) {
    return topMargin + plotHeight - (value / maximum) * plotHeight;
  }

  let markup = "";
  for (let index = 0; index <= 4; index++) {
    const value = (maximum * index) / 4;
    markup += `
      <line x1="${leftMargin}" x2="${width - rightMargin}" y1="${y(value)}" y2="${y(value)}" stroke="#e1e7ef"/>
      <text x="${leftMargin - 8}" y="${y(value) + 4}" text-anchor="end" font-size="11" fill="#647086">${value.toFixed(1)}</text>
    `;
  }
  markup +=
    '<text transform="translate(16,135) rotate(-90)" font-size="12" fill="#647086">−log₁₀(reported p)</text>';

  const pointRadius = data.plot_associations.length > 150 ? 2.5 : 4;
  for (const association of data.plot_associations) {
    const adjusted = association.adjustment_status === "verified_fdr_adjusted";
    const cx = x(association.position);
    const cy = y(association.neg_log10_p);
    const shape = adjusted
      ? `<polygon points="${cx},${cy - 5} ${cx + 5},${cy} ${cx},${cy + 5} ${cx - 5},${cy}" fill="#a64b00"/>`
      : `<circle cx="${cx}" cy="${cy}" r="${pointRadius}" fill="#52616f" fill-opacity="0.75"/>`;
    const label = `${association.trait} · chr${association.chromosome}:${association.position} · reported p=${association.p_value_text} · ${association.statistic_kind} · ${association.review_status} · ${association.relation}`;
    markup += `<g role="img" aria-label="${escapeHtml(label)}">${shape}<title>${escapeHtml(label)}</title></g>`;
  }
  for (let index = 0; index <= 4; index++) {
    const position = Math.round(data.start + ((data.end - data.start) * index) / 4);
    markup += `<text x="${x(position)}" y="210" text-anchor="middle" font-size="11" fill="#647086">${position.toLocaleString()}</text>`;
  }

  const compact = data.neighbors.length > 20 || data.window > 100000;
  const laneEnds = [];
  for (const gene of data.neighbors) {
    const left = x(Math.max(data.start, gene.start));
    const right = x(Math.min(data.end, gene.end));
    let lane = compact ? 0 : findAvailableLane(laneEnds, left);
    if (compact) {
      while (lane < 3 && (laneEnds[lane] ?? -Infinity) >= left - 3) lane++;
    }
    laneEnds[lane] = right;
    const centerY = 250 + lane * (compact ? 22 : 42);
    const color = gene.id === data.gene.id ? "#005a8d" : "#52616f";
    const tip = Math.min(7, (right - left) / 2);
    const points =
      gene.strand === "+"
        ? `${left},${centerY - 5} ${right - tip},${centerY - 5} ${right},${centerY} ${right - tip},${centerY + 5} ${left},${centerY + 5}`
        : `${right},${centerY - 5} ${left + tip},${centerY - 5} ${left},${centerY} ${left + tip},${centerY + 5} ${right},${centerY + 5}`;
    const label = `${gene.id} · ${gene.symbol || "No symbol"} · ${gene.start.toLocaleString()}–${gene.end.toLocaleString()} (${gene.strand}) · ${gene.description || "No description"}`;
    markup += `<g class="clickgene" tabindex="0" role="button" aria-label="Explore ${escapeHtml(label)}" data-gene="${escapeHtml(gene.id)}" data-detail="${escapeHtml(label)}">
        <rect x="${left - 3}" y="${centerY - 9}" width="${Math.max(7, right - left + 6)}" height="18" fill="transparent"/>
        <polygon points="${points}" fill="${color}" stroke="${color}" stroke-width="${gene.id === data.gene.id ? 2 : 0.7}"/>
        ${compact ? "" : `<text x="${Math.min(left, width - rightMargin - 150)}" y="${centerY + 21}" font-size="11" fill="${color}">${escapeHtml(gene.symbol || gene.id)}</text>`}
        <title>${escapeHtml(label)}</title></g>`;
  }
  if (compact) {
    const selected = data.gene;
    const labelX = Math.max(
      leftMargin + 85,
      Math.min(width - rightMargin - 85, x((selected.start + selected.end) / 2)),
    );
    markup += `<line x1="${x((selected.start + selected.end) / 2)}" x2="${labelX}" y1="${250 + 3 * 22 + 12}" y2="345" stroke="#005a8d"/><text x="${labelX}" y="363" text-anchor="middle" font-size="12" fill="#005a8d">${escapeHtml(selected.symbol || selected.id)} · selected gene</text>`;
  }
  svg.setAttribute(
    "viewBox",
    `0 0 ${width} ${compact ? 390 : Math.max(320, 270 + laneEnds.length * 42)}`,
  );
  byId("chart-mode").textContent = compact
    ? `Compact overview · ${data.neighbors.length} genes. Zoom in to see more gene labels.`
    : `${data.neighbors.length} genes · Hover over or focus a gene for details.`;
  byId("gene-detail").textContent = "Select a gene arrow to explore its evidence.";
  byId("zoom-in").disabled = byId("window").selectedIndex === 0;
  byId("zoom-out").disabled =
    byId("window").selectedIndex === byId("window").options.length - 1;

  svg.innerHTML = markup;
  bindGeneControls(svg);
  for (const mark of svg.querySelectorAll("[data-detail]")) {
    mark.onmouseenter = showChartGeneDetail;
    mark.onfocus = showChartGeneDetail;
  }
}

// Open the gene selected through an ambiguity choice, table row or chart arrow.
function chooseGene(event) {
  loadRegion(event.currentTarget.dataset.gene, 1);
}

// Activate a focused chart gene when the user presses Enter or Space.
function chooseGeneWithKeyboard(event) {
  if (event.key === "Enter" || event.key === " ") {
    event.preventDefault();
    chooseGene(event);
  }
}

// Attach mouse and keyboard navigation to gene controls inside one container.
function bindGeneControls(container) {
  for (const control of container.querySelectorAll("[data-gene]")) {
    control.onclick = chooseGene;
    control.onkeydown = chooseGeneWithKeyboard;
  }
}

// Submit a gene search without allowing the browser to reload the page.
function submitSearch(event) {
  event.preventDefault();
  searchGenes();
}

// Reload the selected gene using the current window and trait controls.
function applyFilters() {
  if (selectedGene) loadRegion(selectedGene, 1);
}

// Request the preceding page of regional association records.
function showPreviousPage() {
  loadRegion(selectedGene, currentPage - 1);
}

// Request the following page of regional association records.
function showNextPage() {
  loadRegion(selectedGene, currentPage + 1);
}

// Copy a documented example into the search field and resolve its matching genes.
function chooseExample(event) {
  const { example, exampleStudy } = event.currentTarget.dataset;
  byId("query").value = example;
  if (exampleStudy) {
    if (expressionCatalogLoaded) byId("expression-dataset").value = exampleStudy;
    else initialExpressionDataset = exampleStudy;
  }
  searchGenes();
}

// Restore a bookmarked gene and supported window before starting the first search.
function restoreBookmark() {
  const initial = new URLSearchParams(location.hash.slice(1));
  if (!initial.has("gene")) return;
  byId("query").value = initial.get("gene");
  for (const option of byId("window").options) {
    if (option.value === initial.get("window")) {
      byId("window").value = initial.get("window");
      break;
    }
  }
}

byId("search").onsubmit = submitSearch;
byId("apply").onclick = applyFilters;
byId("reset").onclick = resetFilters;
byId("go-page").onclick = goToPage;
byId("window").onchange = applyFilters;
byId("prev").onclick = showPreviousPage;
byId("next").onclick = showNextPage;
for (const button of document.querySelectorAll("[data-example]")) {
  button.onclick = chooseExample;
}
restoreBookmark();
if (location.pathname !== "/") searchGenes();

// Offer explicit comparison choices for genes in the selected genomic neighborhood.
function renderExpressionControls(neighbors) {
  let html = "<fieldset><legend>Compare neighboring genes (up to 12)</legend>";
  for (const gene of neighbors) {
    html += `<label class="expression-choice"><input type="checkbox" value="${escapeHtml(gene.id)}" ${gene.id === selectedGene ? "checked" : ""}> ${escapeHtml(gene.symbol || gene.id)} <small>${escapeHtml(gene.id)}</small></label>`;
  }
  byId("expression-controls").innerHTML =
    html +
    '</fieldset><button id="compare-expression" type="button">Update expression comparison</button>';
  byId("compare-expression").onclick = loadExpression;
}

// Retrieve baseline expression separately so a dataset problem cannot hide GWAS results.
async function loadExpression() {
  const requestId = ++expressionRequest;
  loadStudyCoverage(requestId);
  if (!expressionCatalogLoaded) {
    try {
      const catalog = await fetchJson("/api/expression/datasets");
      if (requestId !== expressionRequest) return;
      byId("expression-dataset").replaceChildren();
      for (const study of catalog.datasets)
        byId("expression-dataset").add(
          new Option(`${study.citation} · ${study.organ}`, study.id),
        );
      if (initialExpressionDataset) {
        for (const option of byId("expression-dataset").options)
          if (option.value === initialExpressionDataset)
            byId("expression-dataset").value = initialExpressionDataset;
        initialExpressionDataset = "";
      }
      expressionCatalogLoaded = true;
    } catch (error) {
      byId("expression-status").textContent = error.message;
      return;
    }
  }
  const datasetId = byId("expression-dataset").value || "robertson_2026_control";
  const geneIds = [];
  for (const input of byId("expression-controls").querySelectorAll("input:checked"))
    geneIds.push(input.value);
  if (geneIds.includes(selectedGene)) {
    geneIds.splice(geneIds.indexOf(selectedGene), 1);
    geneIds.unshift(selectedGene);
  }
  if (!geneIds.length || geneIds.length > 12) {
    byId("expression-status").textContent = "Select between one and twelve genes.";
    return;
  }
  byId("expression-status").textContent = "Loading baseline expression…";
  byId("expression-content").replaceChildren();
  byId("expression-organ").disabled = true;
  byId("expression-cell").disabled = true;
  expressionData = null;
  try {
    const data = await fetchJson(
      `/api/expression?genes=${encodeURIComponent(geneIds.join(","))}&dataset=${encodeURIComponent(datasetId)}`,
    );
    if (requestId !== expressionRequest) return;
    if (!data.available) {
      byId("expression-status").textContent = data.message;
      return;
    }
    byId("expression-status").textContent =
      `${data.dataset.citation} · ${datasetId === "wang_2025_normal_root" ? "Root only" : data.dataset.organ} · ${data.dataset.condition} · ${data.dataset.n_cells.toLocaleString()} ${data.dataset.cell_unit || "nuclei"} · ${data.dataset.n_samples} replicates`;
    expressionExport = null;
    populateExpressionFilters(data);
    const bookmark = new URLSearchParams(location.hash.slice(1));
    bookmark.set("expression", datasetId);
    history.replaceState(null, "", "#" + bookmark);
    const exportLink = document.createElement("a");
    exportLink.className = "button secondary";
    exportLink.textContent = "Export expression CSV";
    exportLink.href = `/api/expression/export?genes=${encodeURIComponent(geneIds.join(","))}&dataset=${encodeURIComponent(datasetId)}`;
    expressionExport = exportLink;
    byId("expression-content").prepend(exportLink);
    updateToolLinks();
  } catch (error) {
    if (requestId === expressionRequest)
      byId("expression-status").textContent = error.message;
  }
}

// Show source availability separately so supporting tables cannot be mistaken for measured expression.
async function loadStudyCoverage(requestId) {
  const container = byId("study-coverage");
  container.replaceChildren();
  try {
    const data = await fetchJson(
      `/api/expression/sources?gene=${encodeURIComponent(selectedGene)}`,
    );
    if (requestId !== expressionRequest) return;
    if (!data.sources.length) return;
    let html =
      '<details class="study-coverage"><summary>Additional studies and data coverage</summary>';
    for (const source of data.sources) {
      html += `<article><h3><a href="${escapeHtml(source.publication_url)}" target="_blank" rel="noopener noreferrer">${escapeHtml(source.citation)}</a></h3><p><strong>${escapeHtml(source.status)}</strong></p><p>${escapeHtml(source.message)}</p>`;
      for (const mapping of source.gene_mappings || [])
        html += `<p>Author cross-reference: ${escapeHtml(data.gene_id)} ↔ ${escapeHtml(mapping.msu_ids.join(", "))} ↔ ${escapeHtml(mapping.source_gene_id)}.</p>`;
      if (source.mapping_note) html += `<p>${escapeHtml(source.mapping_note)}</p>`;
      html += '<p class="source-links">';
      for (const link of source.links)
        html += `<a href="${escapeHtml(link.url)}" target="_blank" rel="noopener noreferrer">${escapeHtml(link.label)}</a> `;
      html += "</p></article>";
    }
    container.innerHTML = html + "</details>";
  } catch (error) {
    if (requestId === expressionRequest)
      container.textContent =
        "Additional study coverage could not be loaded. Expression results remain available below.";
  }
}

// Render study-specific dot plots and numeric tables without substituting zero for missing genes.
function renderExpression(data) {
  const dataset = data.dataset;
  let html = `<p><a href="${escapeHtml(dataset.publication_url)}" target="_blank" rel="noopener">${escapeHtml(dataset.citation)}</a> · <a href="${escapeHtml(dataset.repository_url)}" target="_blank" rel="noopener">Source data</a></p><details class="reading-guide"><summary>Methods and data coverage</summary><p>${escapeHtml(dataset.metric_description)}</p><p>${escapeHtml(dataset.coverage_note)}</p></details>`;
  if (dataset.annotation_source_url) {
    html += `<p class="notice">Cell annotations: <a href="${escapeHtml(dataset.annotation_source_url)}" target="_blank" rel="noopener">${escapeHtml(dataset.annotation_source_name)}</a> · <a href="${escapeHtml(dataset.annotation_publication_url)}" target="_blank" rel="noopener">Curation publication</a>. ${escapeHtml(dataset.annotation_note)}</p>`;
  }
  const measured = [];
  let resultSummary =
    '<ul class="expression-result-summary" aria-label="Expression results">';
  for (const gene of data.genes) {
    const id = escapeHtml(gene.gene_id);
    if (gene.status !== "measured") {
      resultSummary += `<li><strong>${id}</strong> · Gene not measured in the available data</li>`;
      continue;
    }
    measured.push(gene);
    if (!gene.groups.length) {
      resultSummary += `<li><strong>${id}</strong> · No measurements for the selected cell types</li>`;
      continue;
    }
    let detected = 0;
    for (const group of gene.groups) if (group.n_expressing > 0) detected++;
    const zero = gene.groups.length - detected;
    resultSummary += `<li><strong>${id}</strong> · Detected in ${detected} of ${gene.groups.length} measured cell types${zero ? ` · Zero detected in ${zero}` : ""}</li>`;
  }
  html = resultSummary + "</ul>" + html;
  if (!measured.length) {
    byId("expression-content").innerHTML = html;
    return;
  }
  const first = measured[0];
  const order = [];
  for (const row of first.groups) order.push(row.cell_type);
  for (const group of data.cell_groups)
    if (!order.includes(group.cell_type)) order.push(group.cell_type);
  let maximum = 0;
  for (const gene of measured)
    for (const row of gene.groups) maximum = Math.max(maximum, row.mean_log);
  const width = Math.max(1100, 260 + measured.length * 145);
  const height = 115 + order.length * 32;
  const columnSpacing = (width - 260) / measured.length;
  let svg = `<svg class="expression-plot" viewBox="0 0 ${width} ${height}" role="img" aria-label="Within-study baseline expression; numeric values follow in the table">`;
  for (let col = 0; col < measured.length; col++) {
    const gene = measured[col];
    svg += `<text x="${260 + columnSpacing * (col + 0.5)}" y="22" text-anchor="middle" font-size="11">${escapeHtml(gene.gene_id)}</text>`;
    for (let index = 0; index < order.length; index++) {
      const row = gene.groups.find(findExpressionGroup, order[index]);
      const cx = 260 + columnSpacing * (col + 0.5);
      const cy = 65 + index * 32;
      if (!row) {
        svg += `<text x="${cx}" y="${cy + 4}" text-anchor="middle" font-size="11">NA<title>Feature unavailable for the samples representing this cell type</title></text>`;
        continue;
      }
      const radius = 11 * Math.sqrt(row.percent_expressing / 100);
      const shade = maximum ? row.mean_log / maximum : 0;
      const red = Math.round(225 - shade * 225);
      const green = Math.round(238 - shade * 148);
      const blue = Math.round(247 - shade * 106);
      svg += `<circle cx="${cx}" cy="${cy}" r="${radius}" fill="rgb(${red},${green},${blue})" stroke="#52616f" stroke-width=".7"><title>${escapeHtml(row.cell_type)} · ${escapeHtml(gene.gene_id)} · mean ${row.mean_log.toFixed(3)} · ${row.percent_expressing.toFixed(1)}% detected</title></circle>`;
    }
  }
  for (let index = 0; index < order.length; index++)
    svg += `<text x="225" y="${69 + index * 32}" text-anchor="end" font-size="12">${escapeHtml(order[index])}</text>`;
  svg += "</svg>";
  html += `<p>Larger dots mean a higher percentage of cells expressed the gene. Darker dots mean higher average expression.</p><div class="expression-legend"><svg viewBox="0 0 650 55" role="img" aria-label="Dot area represents percentage detected: ten, fifty and one hundred percent. Color goes from zero to the displayed maximum mean expression."><defs><linearGradient id="expression-gradient"><stop offset="0%" stop-color="#e1eef7"/><stop offset="100%" stop-color="#005a8d"/></linearGradient></defs><circle cx="20" cy="20" r="${11 * Math.sqrt(0.1)}" fill="#52616f"/><text x="35" y="24">10%</text><circle cx="100" cy="20" r="${11 * Math.sqrt(0.5)}" fill="#52616f"/><text x="118" y="24">50%</text><circle cx="190" cy="20" r="11" fill="#52616f"/><text x="208" y="24">100%</text><rect x="300" y="10" width="150" height="20" fill="url(#expression-gradient)"/><text x="295" y="49">0</text><text x="450" y="49" text-anchor="end">${maximum.toFixed(3)}</text><text x="470" y="24">Mean log expression</text></svg></div><div class="chartscroll">${svg}</div>`;
  html += `<details class="reading-guide"><summary>View expression values and ranks</summary><p>Rows are ranked within this study. Interpret small groups with care. In the plot, NA marks unavailable measurements and a blank marks zero detection.</p><div class="tablewrap expression-table"><table><thead><tr><th>Gene</th><th>Cell type</th><th>Tissue system / organ</th><th>Rank</th><th>Mean log expression</th><th>Expressing ${escapeHtml(dataset.cell_unit || "nuclei")}</th><th>Percentage</th><th>Replicates represented</th></tr></thead><tbody>`;
  for (const gene of measured)
    for (const row of gene.groups)
      html += `<tr><td>${escapeHtml(gene.gene_id)}</td><td>${escapeHtml(row.cell_type)}</td><td>${escapeHtml(row.tissue_system)}</td><td>${row.rank}</td><td>${row.mean_log.toFixed(4)}</td><td>${row.n_expressing.toLocaleString()} / ${row.n_cells.toLocaleString()}</td><td>${row.percent_expressing.toFixed(2)}%</td><td>${row.n_samples} / ${dataset.n_samples}</td></tr>`;
  html += "</tbody></table></div></details>";
  byId("expression-content").innerHTML = html;
}

// Match a cell-type summary to the row label requested by the current plot.
function findExpressionGroup(row) {
  return row.cell_type === String(this);
}

byId("expression-dataset").onchange = loadExpression;

// Display the requested tool while preserving one shared gene and browser history.
function showToolPage() {
  const previousPage = document.body.dataset.page;
  const page = location.pathname.slice(1) || "home";
  document.body.dataset.page = page;
  const titles = {
    home: ["RiceLENS", ""],
    variation: ["Genetic Variation", "Find recorded SNP positions near a rice gene."],
    gwas: ["GWAS Associations", "Find published trait associations near a rice gene."],
    expression: ["Cell Type Expression", "See where a rice gene is expressed."],
  };
  const [title, description] = titles[page] || titles.home;
  byId("page-title").textContent = title;
  byId("page-description").textContent = description;
  byId("page-intro").hidden = page === "home";
  document.title =
    page === "home" ? "RiceLENS · Wilkins Lab" : `${title} · RiceLENS · Wilkins Lab`;
  updateToolLinks();
  if (
    selectedGene &&
    previousPage &&
    previousPage !== page &&
    page !== "home" &&
    byId("variant-search").value
  )
    loadRegion(selectedGene, 1);
}

// Carry the current gene, window and expression study into each tool link.
function updateToolLinks() {
  for (const link of document.querySelectorAll("[data-route]")) {
    const target = link.dataset.route;
    const bookmark = new URLSearchParams(location.hash.slice(1));
    if (target === "home" && bookmark.has("gene")) bookmark.set("page", "home");
    else bookmark.delete("page");
    link.href =
      (target === "home" ? "/" : `/${target}`) + (bookmark.size ? "#" + bookmark : "");
    if (target === document.body.dataset.page)
      link.setAttribute("aria-current", "page");
    else link.removeAttribute("aria-current");
  }
}

// Navigate between tools without discarding the search, comparison or filter state.
function navigateTool(event) {
  if (
    event.button !== 0 ||
    event.metaKey ||
    event.ctrlKey ||
    event.shiftKey ||
    event.altKey
  )
    return;
  event.preventDefault();
  history.pushState(null, "", event.currentTarget.href);
  showToolPage();
  if (!selectedGene && document.body.dataset.page !== "home") searchGenes();
  window.scrollTo(0, 0);
  byId("main-content").focus({ preventScroll: true });
}

// Rebuild valid organ and cell-type choices from the selected study response.
function populateExpressionFilters(data) {
  expressionData = data;
  const organs = new Set();
  for (const group of data.cell_groups) organs.add(group.organ || data.dataset.organ);
  byId("expression-organ").replaceChildren(new Option("All available organs", ""));
  for (const organ of organs) byId("expression-organ").add(new Option(organ, organ));
  byId("expression-organ").disabled = organs.size === 0;
  refreshExpressionCells();
}

// Limit cell-type options to the currently selected organ without mixing studies.
function refreshExpressionCells() {
  const organ = byId("expression-organ").value;
  byId("expression-cell").replaceChildren(new Option("All available cell types", ""));
  for (const group of expressionData.cell_groups) {
    if (!organ || (group.organ || expressionData.dataset.organ) === organ)
      byId("expression-cell").add(new Option(group.cell_type, group.cell_type));
  }
  byId("expression-cell").disabled = byId("expression-cell").options.length < 2;
  filterExpression();
}

// Filter displayed contexts while retaining original within-study expression ranks.
function filterExpression() {
  if (!expressionData) return;
  const organ = byId("expression-organ").value;
  const cell = byId("expression-cell").value;
  const filtered = structuredClone(expressionData);
  filtered.cell_groups = [];
  for (const group of expressionData.cell_groups) {
    if (
      (!organ || (group.organ || expressionData.dataset.organ) === organ) &&
      (!cell || group.cell_type === cell)
    )
      filtered.cell_groups.push(group);
  }
  const allowed = new Set();
  for (const group of filtered.cell_groups) allowed.add(group.cell_type);
  for (const gene of filtered.genes) {
    const groups = [];
    for (const group of gene.groups || [])
      if (allowed.has(group.cell_type)) groups.push(group);
    gene.groups = groups;
  }
  renderExpression(filtered);
  if (expressionExport) byId("expression-content").prepend(expressionExport);
}

for (const link of document.querySelectorAll("[data-route]"))
  link.onclick = navigateTool;
window.addEventListener("popstate", showToolPage);
byId("expression-organ").onchange = refreshExpressionCells;
byId("expression-cell").onchange = filterExpression;
showToolPage();

// Display recorded positions and allele fields without inventing variant annotations.
function renderVariant(record) {
  return `<tr><td>chr${record.chromosome}:${Number(record.position).toLocaleString()}</td><td>${escapeHtml(record.relation)} · ${Number(record.distance_bp).toLocaleString()} bp</td><td>${escapeHtml(record.effect_allele || "Not reported")}</td><td>${escapeHtml(record.other_allele || "Not reported")}</td><td>${escapeHtml(record.coordinate_status)}<br><a href="${escapeHtml(record.source_record_url)}" target="_blank" rel="noopener">Source record</a></td></tr>`;
}

// Expose the focused or hovered gene description in a persistent readable detail area.
function showChartGeneDetail(event) {
  byId("gene-detail").textContent = event.currentTarget.dataset.detail;
}

// Narrow the genomic window using the next supported scale around the selected gene.
function zoomInRegion() {
  byId("window").selectedIndex = Math.max(0, byId("window").selectedIndex - 1);
  applyFilters();
}

// Widen the genomic window using the next supported scale around the selected gene.
function zoomOutRegion() {
  byId("window").selectedIndex = Math.min(
    byId("window").options.length - 1,
    byId("window").selectedIndex + 1,
  );
  applyFilters();
}

// Filter the complete neighboring-gene table without removing genes from the plot.
function filterNeighborTable() {
  const query = byId("neighbor-search").value.trim().toLowerCase();
  for (const row of byId("neighbors").rows)
    row.hidden = !row.textContent.toLowerCase().includes(query);
  let count = 0;
  for (const row of byId("neighbors").rows) if (!row.hidden) count++;
  byId("neighbor-search-status").textContent = count
    ? `${count} of ${neighboringGenes.length} genes shown in this window.`
    : "No matching genes in this window. Try a different identifier, symbol or description, or clear the search.";
}

byId("zoom-in").onclick = zoomInRegion;
byId("zoom-out").onclick = zoomOutRegion;
byId("neighbor-search").oninput = suggestNeighborGenes;

// Close stale suggestions without altering the user's typed position or record identifier.
function closeVariantSuggestions() {
  ++variantSuggestionRequest;
  clearTimeout(variantSuggestionTimer);
  byId("variant-suggestions").hidden = true;
  byId("variant-search").setAttribute("aria-expanded", "false");
  byId("variant-search").removeAttribute("aria-activedescendant");
  activeVariantOption = -1;
}

// Debounce all-region suggestions while retaining an explicitly submitted results table.
function queueVariantSuggestions() {
  closeVariantSuggestions();
  const value = byId("variant-search").value.trim();
  byId("variant-search-status").textContent =
    "Choose a suggested record or press Search to apply this text to every matching regional record.";
  if (!value || !selectedGene) return;
  variantSuggestionTimer = setTimeout(fetchVariantSuggestions, 250);
}

// Retrieve suggestions across all matching regional records rather than the displayed page.
async function fetchVariantSuggestions() {
  const requestId = ++variantSuggestionRequest;
  const query = new URLSearchParams(regionQuery);
  query.set("variant_query", byId("variant-search").value);
  try {
    const data = await fetchJson(`/api/variant-suggestions?${query}`);
    if (requestId !== variantSuggestionRequest) return;
    variantOptions = data.suggestions;
    const list = byId("variant-suggestions");
    list.replaceChildren();
    for (let index = 0; index < variantOptions.length; index++) {
      const record = variantOptions[index];
      const option = document.createElement("li");
      option.id = `variant-option-${index}`;
      option.setAttribute("role", "option");
      option.setAttribute("aria-selected", "false");
      option.dataset.index = index;
      option.textContent = `chr${record.chromosome}:${record.position.toLocaleString()} · record ${record.id} · ${record.trait}`;
      option.onclick = chooseVariantSuggestion;
      list.append(option);
    }
    list.hidden = !variantOptions.length;
    byId("variant-search").setAttribute(
      "aria-expanded",
      String(variantOptions.length > 0),
    );
    byId("variant-search-status").textContent = data.total
      ? `${data.total.toLocaleString()} matching records in this region; showing up to 8 suggestions. Use arrow keys and Enter to select, or Search to show all matches.`
      : "No matching records in the selected region and filters. Check the position or record ID, widen the window, or clear the search.";
  } catch (error) {
    if (requestId === variantSuggestionRequest)
      byId("variant-search-status").textContent =
        "Suggestions could not be loaded. You can still submit the search.";
  }
}

// Select the exact association record explicitly chosen from the suggestion list.
function chooseVariantSuggestion(event) {
  selectVariantRecord(Number(event.currentTarget.dataset.index));
}

// Apply a selected record identifier without silently changing unselected user input.
function selectVariantRecord(index) {
  byId("variant-search").value = variantOptions[index].id;
  closeVariantSuggestions();
  byId("variant-search-status").textContent =
    "Showing the selected record across this regional dataset.";
  loadRegion(selectedGene, 1);
}

// Support standard combobox arrow navigation, explicit selection and Escape dismissal.
function navigateVariantSuggestions(event) {
  if (event.key === "Escape" || event.key === "Tab") closeVariantSuggestions();
  const open = !byId("variant-suggestions").hidden;
  if (open && (event.key === "ArrowDown" || event.key === "ArrowUp")) {
    event.preventDefault();
    activeVariantOption =
      (activeVariantOption +
        (event.key === "ArrowDown" ? 1 : -1) +
        variantOptions.length) %
      variantOptions.length;
    for (let index = 0; index < variantOptions.length; index++)
      byId(`variant-option-${index}`).setAttribute(
        "aria-selected",
        String(index === activeVariantOption),
      );
    byId("variant-search").setAttribute(
      "aria-activedescendant",
      `variant-option-${activeVariantOption}`,
    );
  } else if (open && event.key === "Enter" && activeVariantOption >= 0) {
    event.preventDefault();
    selectVariantRecord(activeVariantOption);
  }
}

// Apply manually typed variant text to the complete region and reset pagination.
function submitVariantSearch(event) {
  event.preventDefault();
  closeVariantSuggestions();
  byId("variant-search-status").textContent =
    "Search applied to all matching records in the selected region. Clear to restore regional results.";
  if (selectedGene) loadRegion(selectedGene, 1);
}

// Restore all regional records without changing the selected gene or scientific filters.
function clearVariantSearch() {
  byId("variant-search").value = "";
  closeVariantSuggestions();
  byId("variant-search-status").textContent =
    "Search cleared. Showing all records matching the region and other filters.";
  if (selectedGene) loadRegion(selectedGene, 1);
}

byId("variant-search").oninput = queueVariantSuggestions;
byId("variant-search").onkeydown = navigateVariantSuggestions;
byId("variant-search-form").onsubmit = submitVariantSearch;
byId("clear-variant-search").onclick = clearVariantSearch;

// Dismiss the suggestion dropdown when the user clicks outside its search form.
function dismissVariantSuggestions(event) {
  if (!byId("variant-search-form").contains(event.target)) closeVariantSuggestions();
}

document.addEventListener("click", dismissVariantSuggestions);

// Close neighboring-gene suggestions without changing the entered query or table.
function closeNeighborSuggestions() {
  byId("neighbor-suggestions").hidden = true;
  byId("neighbor-search").setAttribute("aria-expanded", "false");
  byId("neighbor-search").removeAttribute("aria-activedescendant");
  activeNeighborOption = -1;
}

// Suggest matching IDs, symbols and descriptions from every gene in the current window.
function suggestNeighborGenes() {
  closeNeighborSuggestions();
  const query = byId("neighbor-search").value.trim().toLowerCase();
  const list = byId("neighbor-suggestions");
  list.replaceChildren();
  neighborOptions = [];
  if (!query) return;
  let count = 0;
  for (const gene of neighboringGenes) {
    if (
      `${gene.id} ${gene.symbol || ""} ${gene.description || ""}`
        .toLowerCase()
        .includes(query)
    ) {
      count++;
      if (neighborOptions.length < 8) neighborOptions.push(gene);
    }
  }
  for (let index = 0; index < neighborOptions.length; index++) {
    const gene = neighborOptions[index];
    const option = document.createElement("li");
    option.id = `neighbor-option-${index}`;
    option.setAttribute("role", "option");
    option.setAttribute("aria-selected", "false");
    option.dataset.index = index;
    option.textContent = `${gene.id} · ${gene.symbol || "No symbol"} · ${gene.description || "No description"}`;
    option.onclick = chooseNeighborSuggestion;
    list.append(option);
  }
  list.hidden = neighborOptions.length === 0;
  byId("neighbor-search").setAttribute(
    "aria-expanded",
    String(neighborOptions.length > 0),
  );
  byId("neighbor-search-status").textContent = count
    ? `${count} matching genes; showing up to 8 suggestions. Choose one or Search to show all matches.`
    : "No matching genes in this window. Try another ID, symbol or description, or widen the window.";
}

// Filter the table to the exact gene explicitly chosen from the suggestions.
function selectNeighborSuggestion(index) {
  byId("neighbor-search").value = neighborOptions[index].id;
  closeNeighborSuggestions();
  filterNeighborTable();
}

// Apply the neighboring gene selected by clicking a suggestion.
function chooseNeighborSuggestion(event) {
  selectNeighborSuggestion(Number(event.currentTarget.dataset.index));
}

// Support arrow keys, Enter selection and Escape or Tab dismissal for gene suggestions.
function navigateNeighborSuggestions(event) {
  if (event.key === "Escape" || event.key === "Tab") closeNeighborSuggestions();
  const open = !byId("neighbor-suggestions").hidden;
  if (open && (event.key === "ArrowDown" || event.key === "ArrowUp")) {
    event.preventDefault();
    activeNeighborOption =
      (activeNeighborOption +
        (event.key === "ArrowDown" ? 1 : -1) +
        neighborOptions.length) %
      neighborOptions.length;
    for (let index = 0; index < neighborOptions.length; index++)
      byId(`neighbor-option-${index}`).setAttribute(
        "aria-selected",
        String(index === activeNeighborOption),
      );
    byId("neighbor-search").setAttribute(
      "aria-activedescendant",
      `neighbor-option-${activeNeighborOption}`,
    );
  } else if (open && event.key === "Enter" && activeNeighborOption >= 0) {
    event.preventDefault();
    selectNeighborSuggestion(activeNeighborOption);
  }
}

// Apply the typed neighboring-gene query only when the user submits the search.
function submitNeighborSearch(event) {
  event.preventDefault();
  closeNeighborSuggestions();
  filterNeighborTable();
}

// Restore the full neighboring-gene table without changing the selected genomic region.
function clearNeighborSearch() {
  byId("neighbor-search").value = "";
  closeNeighborSuggestions();
  filterNeighborTable();
}

// Dismiss neighboring-gene suggestions when focus moves outside through a mouse click.
function dismissNeighborSuggestions(event) {
  if (!byId("neighbor-search-form").contains(event.target)) closeNeighborSuggestions();
}

byId("neighbor-search").onkeydown = navigateNeighborSuggestions;
byId("neighbor-search-form").onsubmit = submitNeighborSearch;
byId("clear-neighbor-search").onclick = clearNeighborSearch;
document.addEventListener("click", dismissNeighborSuggestions);
