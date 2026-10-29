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
      `${data.matches.length} matching genes${limitNotice}. Select one to continue.`,
    );
    for (const gene of data.matches) {
      const button = document.createElement("button");
      button.textContent = `${gene.id} · ${gene.symbol || "unnamed"} · chr${gene.chromosome}`;
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
  const allele = association.effect_allele
    ? `Effect allele: ${escapeHtml(association.effect_allele)}`
    : "Effect allele not reported";
  const frequency = association.effect_allele_frequency
    ? `Frequency: ${escapeHtml(association.effect_allele_frequency)}`
    : "Frequency not reported";
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
      <td>${escapeHtml(association.trait)}<small>${escapeHtml(association.ontology_id)}</small></td>
      <td>
        chr${association.chromosome}:${association.position.toLocaleString()}
        <small>Source row ${association.source_row} · record ${escapeHtml(association.id)}</small>
      </td>
      <td>
        ${escapeHtml(association.relation)}
        <small>${association.distance_bp.toLocaleString()} bp from gene span</small>
      </td>
      <td>${escapeHtml(association.p_value_text)}</td>
      <td>${allele}<small>${frequency}<br>Effect size: not supplied in export${variance}</small></td>
      <td>
        ${escapeHtml(association.source_study_id)} · ${escapeHtml(association.model)}
        <small>
          N = ${escapeHtml(association.sample_size)} · ${escapeHtml(association.environment)}
          · ${escapeHtml(association.tissue)}
        </small>
        ${publication}
        <details>
          <summary>Coordinate provenance</summary>
          Atlas reference: ${escapeHtml(association.source_reference)}<br>
          Original reported location: ${escapeHtml(association.reported_location)}<br>
          <a href="${escapeHtml(association.source_url)}" target="_blank" rel="noopener">Source dataset</a>
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
    const query = new URLSearchParams({
      gene: geneId,
      window: byId("window").value,
      trait: byId("trait").value,
      page,
    });
    const data = await fetchJson(`/api/region?${query}`);
    if (requestId !== latestRequest) return;

    selectedGene = geneId;
    currentPage = page;
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

    const traitNotice = data.trait_filter ? ` matching “${data.trait_filter}”` : "";
    byId("counts").textContent =
      `${totalAssociations.toLocaleString()} curated association records${traitNotice}. ` +
      "Records may share SNP positions, traits or study samples; counts are not independent replications.";
    byId("associations").innerHTML = data.associations.length
      ? data.associations.map(renderAssociation).join("")
      : '<tr><td colspan="6">No records in the imported dataset for this window and trait filter. This does not establish absence of association.</td></tr>';
    byId("neighbors").innerHTML = data.neighbors.map(renderNeighbor).join("");
    bindGeneControls(byId("neighbors"));

    byId("prev").disabled = currentPage === 1;
    byId("next").disabled = currentPage * 100 >= totalAssociations;
    byId("page").textContent =
      `Page ${currentPage} of ${Math.max(1, Math.ceil(totalAssociations / 100))}`;
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
      <line x1="${leftMargin}" x2="${width - rightMargin}" y1="${y(value)}" y2="${y(value)}" stroke="#dce3d9"/>
      <text x="${leftMargin - 8}" y="${y(value) + 4}" text-anchor="end" font-size="11" fill="#586c63">${value.toFixed(1)}</text>
    `;
  }
  markup +=
    '<text transform="translate(16,135) rotate(-90)" font-size="12" fill="#586c63">−log₁₀(p)</text>';

  for (const association of data.plot_associations) {
    const color = association.relation === "overlaps gene" ? "#c9782b" : "#387a65";
    markup += `
      <circle cx="${x(association.position)}" cy="${y(association.neg_log10_p)}" r="4" fill="${color}" opacity=".65">
        <title>${escapeHtml(association.trait)} · chr${association.chromosome}:${association.position} · p=${escapeHtml(association.p_value_text)}</title>
      </circle>
    `;
  }
  for (let index = 0; index <= 4; index++) {
    const position = Math.round(data.start + ((data.end - data.start) * index) / 4);
    markup += `<text x="${x(position)}" y="210" text-anchor="middle" font-size="11" fill="#586c63">${position.toLocaleString()}</text>`;
  }

  const laneEnds = [];
  for (const gene of data.neighbors) {
    const left = x(Math.max(data.start, gene.start));
    const right = x(Math.min(data.end, gene.end));
    const lane = findAvailableLane(laneEnds, left);
    laneEnds[lane] = right;
    const centerY = 250 + lane * 42;
    const color = gene.id === data.gene.id ? "#c9782b" : "#387a65";
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
byId("window").onchange = applyFilters;
byId("prev").onclick = showPreviousPage;
byId("next").onclick = showNextPage;
for (const button of document.querySelectorAll("[data-example]")) {
  button.onclick = chooseExample;
}
restoreBookmark();
searchGenes();
