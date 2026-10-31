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

let selectedGene = "";
let currentPage = 1;
let totalAssociations = 0;
let latestRequest = 0;

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
  const requestId = ++latestRequest;
  setStatus("Loading gene and regional evidence…");

  try {
    if (geneId !== selectedGene) {
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
    });
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

    const plotNotice =
      totalAssociations > data.plot_count
        ? " strongest records; table and CSV include all matching records"
        : "";
    byId("region").textContent =
      `Chr ${data.gene.chromosome}: ${data.start.toLocaleString()}–${data.end.toLocaleString()} · ` +
      `${data.neighbors.length} genes · ${totalAssociations.toLocaleString()} association records. ` +
      `Plot shows ${data.plot_count.toLocaleString()}${plotNotice}.`;

    const traitNotice = data.trait_filter
      ? ` matching “${data.trait_filter}”`
      : "";
    byId("counts").textContent =
      `${totalAssociations.toLocaleString()} evidence records${traitNotice} from ${data.source_record_total.toLocaleString()} matching source rows; ${data.repeated_rows_hidden} exact repeats grouped. ` +
      "Records may share SNP positions, traits or study samples; counts are not independent replications.";
    byId("associations").innerHTML = data.associations.length
      ? data.associations.map(renderAssociation).join("")
      : '<tr><td colspan="6">No records in the imported dataset for this window and trait filter. This does not establish absence of association.</td></tr>';
    byId("neighbors").innerHTML = data.neighbors.map(renderNeighbor).join("");
    bindGeneControls(byId("neighbors"));

    byId("prev").disabled = currentPage === 1;
    byId("next").disabled = currentPage * data.page_size >= totalAssociations;
    byId("page").textContent =
      `Page ${currentPage} of ${Math.max(1, Math.ceil(totalAssociations / data.page_size))}`;
    byId("page-number").value = currentPage;
    byId("page-number").max = data.page_count;
    byId("export").href = `/api/export?${query}`;
    drawNeighborhood(data);
    history.replaceState(
      null,
      "",
      "#" + new URLSearchParams({ gene: geneId, window: data.window }),
    );
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
    byId("study").add(
      new Option(`${item.source_study_id} · ${item.model}`, item.id),
    );
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

  for (const association of data.plot_associations) {
    const adjusted = association.adjustment_status === "verified_fdr_adjusted";
    const cx = x(association.position);
    const cy = y(association.neg_log10_p);
    const shape = adjusted
      ? `<polygon points="${cx},${cy - 5} ${cx + 5},${cy} ${cx},${cy + 5} ${cx - 5},${cy}" fill="#a64b00"/>`
      : `<circle cx="${cx}" cy="${cy}" r="4" fill="#52616f"/>`;
    const label = `${association.trait} · chr${association.chromosome}:${association.position} · reported p=${association.p_value_text} · ${association.statistic_kind} · ${association.review_status} · ${association.relation}`;
    markup += `<g tabindex="0" role="img" aria-label="${escapeHtml(label)}">${shape}<title>${escapeHtml(label)}</title></g>`;
  }
  for (let index = 0; index <= 4; index++) {
    const position = Math.round(
      data.start + ((data.end - data.start) * index) / 4,
    );
    markup += `<text x="${x(position)}" y="210" text-anchor="middle" font-size="11" fill="#647086">${position.toLocaleString()}</text>`;
  }

  const laneEnds = [];
  for (const gene of data.neighbors) {
    const left = x(Math.max(data.start, gene.start));
    const right = x(Math.min(data.end, gene.end));
    const lane = findAvailableLane(laneEnds, left);
    laneEnds[lane] = right;
    const centerY = 250 + lane * 42;
    const color = gene.id === data.gene.id ? "#005a8d" : "#52616f";
    const tip = Math.min(7, (right - left) / 2);
    const points =
      gene.strand === "+"
        ? `${left},${centerY - 5} ${right - tip},${centerY - 5} ${right},${centerY} ${right - tip},${centerY + 5} ${left},${centerY + 5}`
        : `${right},${centerY - 5} ${left + tip},${centerY - 5} ${left},${centerY} ${left + tip},${centerY + 5} ${right},${centerY + 5}`;
    markup += `
      <g class="clickgene" tabindex="0" role="button" aria-label="Explore ${escapeHtml(gene.id)}" data-gene="${escapeHtml(gene.id)}">
        <polygon points="${points}" fill="${color}"/>
        <text x="${left}" y="${centerY + 21}" font-size="11" fill="${color}">${escapeHtml(gene.symbol || gene.id)}</text>
        <title>${escapeHtml(gene.id)} · ${gene.start}–${gene.end} (${gene.strand})</title>
      </g>
    `;
  }
  svg.setAttribute(
    "viewBox",
    `0 0 ${width} ${Math.max(320, 270 + laneEnds.length * 42)}`,
  );
  svg.innerHTML = markup;
  bindGeneControls(svg);
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
  byId("query").value = event.currentTarget.dataset.example;
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
searchGenes();
