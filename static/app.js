// Find a page element by its stable identifier.
function byId(id) {
  return document.getElementById(id);
}
// Escape source text before rendering any HTML or SVG content.
function escapeHtml(value) {
  return String(value ?? "").replace(/[&<>"']/g, escapeCharacter);
}
// Replace an HTML-sensitive character with its safe entity.
function escapeCharacter(character) {
  return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[
    character
  ];
}
const routes = new Set([
  "home",
  "locus",
  "gwas",
  "expression",
  "homolog",
  "download",
  "tutorial",
  "contact",
]);
let selectedGene = "";
let currentPage = 1;
let sortOrder = "p_asc";
let latestRequest = 0;
let suggestionRequest = 0;
let searchRequest = 0;
let suggestionTimer;
let suggestionOptions = [];
let activeSuggestion = -1;
let neighboringGenes = [];
let regionData = null;
let loading = false;
// Read JSON and turn API errors into a human-readable message.
async function fetchJson(path) {
  const response = await fetch(path);
  const data = await response.json();
  if (!response.ok)
    throw new Error(data.error || "The request could not be completed.");
  return data;
}
// Announce search progress and errors without shifting keyboard focus.
function setStatus(text, error = false) {
  byId("message").textContent = text;
  byId("message").className = error ? "error" : "";
  byId("evidence-message").textContent = error ? text : "";
}
// Return the canonical tool name while keeping old variation bookmarks usable.
function currentTool() {
  const name = location.pathname.replace(/^\/+|\/+$/g, "") || "home";
  return name === "variation" ? "locus" : routes.has(name) ? name : "home";
}
// Serialize the complete search context for reloads, bookmarks and exports.
function searchParameters(page = currentPage) {
  return new URLSearchParams({
    gene: selectedGene,
    window: byId("window").value,
    trait: byId("trait").value,
    maximum_p: byId("maximum-p").value,
    position: byId("position").value,
    sort: sortOrder,
    page,
    include_repeats: String(byId("include-repeats").checked),
  });
}
// Update navigation and downloads with one consistent selected region and filter state.
function updateToolLinks() {
  const parameters = selectedGene ? searchParameters() : new URLSearchParams();
  const hash = parameters.size ? "#" + parameters : "";
  for (const link of document.querySelectorAll("[data-route]")) {
    const target = link.dataset.route;
    link.href = (target === "home" ? "/" : "/" + target) + hash;
    if (target === currentTool()) link.setAttribute("aria-current", "page");
    else link.removeAttribute("aria-current");
  }
  byId("download-empty").hidden = Boolean(selectedGene);
  byId("download-ready").hidden = !selectedGene || loading || !regionData;
  if (selectedGene && regionData) {
    const caption = `chr${regionData.gene.chromosome}:${regionData.start.toLocaleString()}–${regionData.end.toLocaleString()} · ±${Number(byId("window").value).toLocaleString()} bp`;
    byId("download-context").textContent =
      `Current selection: ${selectedGene} · ${caption}`;
    byId("download-genes").href =
      "/api/genes/export?" +
      new URLSearchParams({ gene: selectedGene, window: byId("window").value });
    byId("download-gwas").href = "/api/export?" + parameters;
    byId("export").href = "/api/export?" + parameters;
  }
}
// Display one page while retaining the same sidebar DOM and search controls.
function showToolPage() {
  const tool = currentTool();
  document.body.dataset.page = tool;
  for (const panel of document.querySelectorAll("[data-panel]"))
    panel.hidden = panel.dataset.panel !== tool;
  byId("search-sidebar").hidden = !["locus", "gwas", "expression"].includes(tool);
  document
    .querySelector(".workspace")
    .classList.toggle("with-sidebar", !byId("search-sidebar").hidden);
  const titles = {
    home: "Home",
    locus: "Locus and Genes",
    gwas: "GWAS",
    expression: "Cell Expression",
    homolog: "Homolog Search",
    download: "Download",
    tutorial: "Tutorial",
    contact: "Contact Us",
  };
  document.title = `${titles[tool]} · RICELens · Wilkins Lab`;
  closeSuggestions();
  updateToolLinks();
  if (tool === "gwas" && regionData) drawNeighborhood(regionData);
  if (tool === "expression") loadWangAtlas();
}
// Navigate locally without fetching a second copy of the shared application state.
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
  window.scrollTo(0, 0);
  byId("main-content").focus({ preventScroll: true });
}
// Close suggestions and invalidate any response from an earlier typed query.
function closeSuggestions() {
  ++suggestionRequest;
  clearTimeout(suggestionTimer);
  byId("search-suggestions").hidden = true;
  byId("query").setAttribute("aria-expanded", "false");
  byId("query").removeAttribute("aria-activedescendant");
  suggestionOptions = [];
  activeSuggestion = -1;
}
// Search gene names or validated positions without automatically resolving ambiguous aliases.
async function searchGenes(event) {
  if (event) event.preventDefault();
  closeSuggestions();
  const query = byId("query").value.trim();
  if (!query) {
    setStatus("Enter a gene name, ID or genomic position.", true);
    return;
  }
  const requestId = ++searchRequest;
  setStatus("Searching…");
  byId("matches").replaceChildren();
  try {
    const data = await fetchJson("/api/search?" + new URLSearchParams({ q: query }));
    if (requestId !== searchRequest) return;
    if (data.match_type === "locus") return loadRegion(data.selector, 1);
    if (data.matches.length === 1) return loadRegion(data.matches[0].id, 1);
    if (!data.matches.length) {
      setStatus(
        "No matching gene. Try a RAP/MSU identifier or a genomic position.",
        true,
      );
      return;
    }
    setStatus(
      `${data.matches.length}${data.truncated ? "+" : ""} matching genes. Compare the descriptions and select one.`,
    );
    for (const gene of data.matches) {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "gene-choice";
      button.dataset.gene = gene.id;
      button.innerHTML = `<strong>${escapeHtml(gene.symbol || gene.id)}</strong><span>${escapeHtml(gene.id)} · chr${gene.chromosome}</span><span>${escapeHtml(gene.description)}</span>`;
      button.onclick = chooseGene;
      byId("matches").append(button);
    }
  } catch (error) {
    if (requestId === searchRequest) setStatus(error.message, true);
  }
}
// Request short gene suggestions after typing pauses while retaining exact identifiers.
function queueSuggestions() {
  closeSuggestions();
  ++searchRequest;
  const query = byId("query").value.trim();
  if (query.length < 2 || query.includes(":")) return;
  suggestionTimer = setTimeout(loadSuggestions, 220);
}
// Render up to six identifier suggestions with descriptions for distinguishing similar names.
async function loadSuggestions() {
  const query = byId("query").value.trim();
  const request = ++suggestionRequest;
  try {
    const data = await fetchJson("/api/search?" + new URLSearchParams({ q: query }));
    if (request !== suggestionRequest) return;
    suggestionOptions = data.matches.slice(0, 6);
    byId("search-suggestions").innerHTML = suggestionOptions
      .map(renderSuggestion)
      .join("");
    byId("search-suggestions").hidden = !suggestionOptions.length;
    byId("query").setAttribute(
      "aria-expanded",
      String(Boolean(suggestionOptions.length)),
    );
    for (const option of byId("search-suggestions").querySelectorAll("[data-index]"))
      option.onclick = chooseSuggestion;
  } catch (error) {
    if (request === suggestionRequest) closeSuggestions();
  }
}
// Format an individual suggestion as an accessible listbox option.
function renderSuggestion(gene, index) {
  return `<li id="suggestion-${index}" role="option" aria-selected="false" data-index="${index}"><strong>${escapeHtml(gene.symbol || gene.id)}</strong><span>${escapeHtml(gene.id)}</span><small>${escapeHtml(gene.description)}</small></li>`;
}
// Apply an explicitly clicked identifier suggestion.
function chooseSuggestion(event) {
  selectSuggestion(Number(event.currentTarget.dataset.index));
}
// Load a suggested gene only after an explicit click or keyboard confirmation.
function selectSuggestion(index) {
  const gene = suggestionOptions[index];
  if (!gene) return;
  byId("query").value = gene.id;
  closeSuggestions();
  loadRegion(gene.id, 1);
}
// Support arrow navigation, Enter selection and Escape dismissal in the suggestion list.
function navigateSuggestions(event) {
  if (event.key === "Escape") {
    closeSuggestions();
    return;
  }
  if (!suggestionOptions.length) return;
  if (event.key === "ArrowDown" || event.key === "ArrowUp") {
    event.preventDefault();
    activeSuggestion =
      (activeSuggestion +
        (event.key === "ArrowDown" ? 1 : -1) +
        suggestionOptions.length) %
      suggestionOptions.length;
    for (const option of byId("search-suggestions").children)
      option.setAttribute(
        "aria-selected",
        String(Number(option.dataset.index) === activeSuggestion),
      );
    byId("query").setAttribute(
      "aria-activedescendant",
      `suggestion-${activeSuggestion}`,
    );
    byId(`suggestion-${activeSuggestion}`).scrollIntoView({ block: "nearest" });
  } else if (event.key === "Enter" && activeSuggestion >= 0) {
    event.preventDefault();
    selectSuggestion(activeSuggestion);
  } else if (event.key === "Tab") closeSuggestions();
}
// Display selected reference coordinates without implying a coordinate query is a gene.
function renderSelection(data) {
  const gene = data.gene;
  byId("selection").hidden = false;
  byId("selected-context").innerHTML =
    `<strong>${escapeHtml(gene.id || data.selector)}</strong>${gene.id ? `<p>${escapeHtml(gene.symbol || "No display name")}</p>` : ""}<p>chr${gene.chromosome}:${gene.start.toLocaleString()}–${gene.end.toLocaleString()}${gene.id ? ` · ${gene.strand} strand` : ""}</p>`;
  byId("query").value = data.selector;
  byId("position").options[2].textContent =
    data.selection_type === "gene" ? "Upstream of gene" : "Lower genomic coordinate";
  byId("position").options[3].textContent =
    data.selection_type === "gene" ? "Downstream of gene" : "Higher genomic coordinate";
}
// Display source missing-value markers clearly while preserving raw fields in the API and exports.
function sourceText(value) {
  return value == null ||
    ["", "na", "n/a", "null", "none"].includes(String(value).trim().toLowerCase())
    ? "Not reported"
    : String(value);
}
// Render concise source evidence with optional details instead of crowding every table cell.
function renderAssociation(a) {
  const publication = a.publication_url
    ? `<a href="${escapeHtml(a.publication_url)}" target="_blank" rel="noopener noreferrer">${escapeHtml(a.pmid ? "PMID " + a.pmid : a.title)}</a>`
    : "Not reported";
  const note = a.curation_note
    ? `<p class="curation-warning">${escapeHtml(a.curation_note)}</p>`
    : "";
  return `<tr><td>${escapeHtml(a.trait)}<small>${escapeHtml(a.evidence_type)}${a.qtl_label ? ` · QTL: ${escapeHtml(a.qtl_label)}` : ""}</small></td><td class="coordinates">chr${a.chromosome}:${a.position.toLocaleString()}<small>${escapeHtml(a.relation)} · ${a.distance_bp.toLocaleString()} bp</small></td><td class="p-value">${escapeHtml(a.p_value_text)}<small>${a.adjustment_status === "verified_fdr_adjusted" ? "◆ Verified FDR-adjusted" : "● Adjustment unverified"}</small></td><td>${publication}<small>Study ${escapeHtml(a.source_study_id || "not reported")}</small><details><summary>Source details</summary><p>${escapeHtml(a.title)}</p><p>${escapeHtml(a.coordinate_status)}</p><p>Original location: ${escapeHtml(a.reported_location || "Not reported")}</p><p>${escapeHtml(a.review_status)} · ${a.source_record_count} source row(s)</p>${note}<a href="${escapeHtml(a.source_record_url)}" target="_blank" rel="noopener noreferrer">Original source record</a></details></td><td>${escapeHtml(sourceText(a.model))}</td><td>${escapeHtml(sourceText(a.sample_size))}</td><td>${escapeHtml(sourceText(a.environment))}<small>${escapeHtml(sourceText(a.experimental_condition))}</small></td></tr>`;
}
// Show aliases as expandable source information without assigning shared symbols to a single locus.
function renderLocusGene(gene) {
  const variants = gene.variants;
  return `<tr class="${gene.id === selectedGene ? "selected-row" : ""}"><td><button class="text-button gene-link" data-gene="${escapeHtml(gene.id)}">${escapeHtml(gene.id)}</button><small>${escapeHtml(gene.symbol || "No display name")}${gene.id === selectedGene ? " · selected" : ""}</small><details><summary>Description</summary><p>${escapeHtml(gene.description)}</p></details></td><td>${gene.chromosome}</td><td>${gene.start.toLocaleString()}</td><td>${gene.end.toLocaleString()}</td><td><details><summary>${gene.aliases.length} identifiers</summary>${gene.aliases.map(escapeHtml).join("<br>")}</details></td><td>${variants.length ? `<details><summary>${new Set(variants.map(variantPosition)).size} SNP positions</summary>${variants.map(renderLocusVariant).join("")}</details>` : "None recorded"}</td><td>Not reported</td></tr>`;
}
// Extract a position for counting distinct source-observed sites.
function variantPosition(variant) {
  return variant.position;
}
// Preserve effect-allele terminology when verified reference and alternate alleles are absent.
function renderLocusVariant(variant) {
  return `<p>${variant.position.toLocaleString()}${variant.effect_allele ? `<small>Reported effect allele: ${escapeHtml(variant.effect_allele)}</small>` : ""}</p>`;
}
// Render a neighbouring gene as an explicit navigation action.
function renderNeighbor(gene) {
  return `<tr class="${gene.id === selectedGene ? "selected-row" : ""}"><td><button class="text-button gene-link" data-gene="${escapeHtml(gene.id)}">${escapeHtml(gene.id)}</button>${gene.id === selectedGene ? "<small>Selected gene</small>" : ""}</td><td>${escapeHtml(gene.symbol || "—")}</td><td>chr${gene.chromosome}:${gene.start.toLocaleString()}–${gene.end.toLocaleString()} (${gene.strand})</td><td>${escapeHtml(gene.description)}</td></tr>`;
}
// Filter the complete regional gene list independently from the association filters.
function filterNeighborTable() {
  const query = byId("neighbor-search").value.trim().toLowerCase();
  const genes = neighboringGenes.filter(matchesNeighbor);
  // Compare one gene with the current literal name, identifier or description query.
  function matchesNeighbor(gene) {
    return `${gene.id} ${gene.symbol || ""} ${gene.description || ""}`
      .toLowerCase()
      .includes(query);
  }
  byId("neighbors").innerHTML = genes.length
    ? genes.map(renderNeighbor).join("")
    : '<tr><td colspan="4">No matching gene in this window.</td></tr>';
  byId("neighbor-count").textContent =
    `${genes.length} of ${neighboringGenes.length} genes`;
  bindGeneControls(byId("neighbors"));
}
// Present summary statistics computed across all filtered records rather than the visible page.
function renderSummary(data) {
  const summary = data.summary;
  const top = summary.top_traits.map(renderTopTrait).join("");
  byId("summary").innerHTML =
    `<strong class="summary-count">${summary.records.toLocaleString()}</strong><p>GWAS records in this window${byId("trait").value ? " matching your filters" : ""}.</p><hr><h3>Most represented trait${summary.top_traits.length > 1 ? "s (tied)" : ""}</h3>${top || "<p>No records</p>"}<h3>Smallest reported P-value</h3><p class="p-value">${summary.smallest_p === null ? "No records" : summary.smallest_p.toExponential(3)}</p><p class="help">Record counts do not imply independent replication.</p>`;
}
// Render a leading trait and its exact record count, including ties.
function renderTopTrait(trait) {
  return `<p>${escapeHtml(trait.trait)}<small>${trait.count.toLocaleString()} records</small></p>`;
}
// Fetch the full regional gene table and filtered association page as one consistent selection.
async function loadRegion(selector, page = 1) {
  closeSuggestions();
  const request = ++latestRequest;
  setStatus("Loading regional evidence…");
  loading = true;
  byId("main-content").setAttribute("aria-busy", "true");
  byId("download-ready").hidden = true;
  const parameters = searchParameters(page);
  parameters.set("gene", selector);
  try {
    const [data, genes] = await Promise.all([
      fetchJson("/api/region?" + parameters),
      fetchJson(
        "/api/genes?" +
          new URLSearchParams({ gene: selector, window: byId("window").value }),
      ),
    ]);
    if (request !== latestRequest) return;
    selectedGene = data.selector;
    currentPage = data.page;
    regionData = data;
    loading = false;
    neighboringGenes = data.neighbors;
    renderSelection(data);
    const caption = `chr${data.gene.chromosome}:${data.start.toLocaleString()}–${data.end.toLocaleString()} · IRGSP-1.0`;
    byId("locus-region").textContent = caption;
    byId("gene-count").textContent = `${genes.genes.length} genes`;
    byId("genes").innerHTML = genes.genes.length
      ? genes.genes.map(renderLocusGene).join("")
      : '<tr><td colspan="7">No annotated genes overlap this interval.</td></tr>';
    bindGeneControls(byId("genes"));
    byId("region").textContent =
      `${caption} · ${data.neighbors.length} genes · Plot shows ${data.plot_count} strongest of ${data.total_associations.toLocaleString()} matching records.`;
    byId("trait-options").replaceChildren();
    for (const trait of data.trait_options)
      byId("trait-options").append(new Option(trait, trait));
    byId("counts").textContent =
      `${data.total_associations.toLocaleString()} evidence records · ${data.source_record_total.toLocaleString()} source rows · ${data.repeated_rows_hidden.toLocaleString()} exact repeats grouped. Showing ${data.associations.length} per page.`;
    byId("associations").innerHTML = data.associations.length
      ? data.associations.map(renderAssociation).join("")
      : '<tr><td colspan="7">No matching published records. Try a wider window or reset the filters; this does not establish absence of association.</td></tr>';
    byId("prev").disabled = data.page <= 1;
    byId("next").disabled = data.page >= data.page_count;
    byId("page").textContent = `Page ${data.page} of ${data.page_count}`;
    byId("page-number").value = data.page;
    byId("page-number").max = data.page_count;
    byId("p-heading").setAttribute(
      "aria-sort",
      sortOrder === "p_asc" ? "ascending" : "descending",
    );
    byId("sort-p").textContent = sortOrder === "p_asc" ? "P-value ↑" : "P-value ↓";
    byId("matches").replaceChildren();
    for (const panel of document.querySelectorAll(".requires-selection"))
      panel.hidden = false;
    for (const panel of document.querySelectorAll(".empty-selection"))
      panel.hidden = true;
    renderSummary(data);
    if (currentTool() === "expression") loadWangAtlas();
    filterNeighborTable();
    drawNeighborhood(data);
    setStatus("");
    history.replaceState(null, "", location.pathname + "#" + searchParameters());
    updateToolLinks();
  } catch (error) {
    if (request === latestRequest) {
      loading = false;
      setStatus(
        error.message +
          (regionData
            ? " Previous results remain unchanged until a search succeeds."
            : ""),
        true,
      );
      byId("download-ready").hidden = !regionData;
    }
  } finally {
    if (request === latestRequest) byId("main-content").removeAttribute("aria-busy");
  }
}
// Re-centre every regional view on the clicked gene while retaining the chosen window and filters.
function chooseGene(event) {
  loadRegion(event.currentTarget.dataset.gene, 1);
}
// Activate an SVG gene control with standard keyboard actions.
function chooseGeneWithKeyboard(event) {
  if (event.key === "Enter" || event.key === " ") {
    event.preventDefault();
    chooseGene(event);
  }
}
// Bind explicit gene-selection actions to table buttons and chart arrows.
function bindGeneControls(container) {
  for (const control of container.querySelectorAll("[data-gene]")) {
    control.onclick = chooseGene;
    if (control.tagName.toLowerCase() !== "button")
      control.onkeydown = chooseGeneWithKeyboard;
  }
}
// Announce an individual chart gene when it receives pointer or keyboard focus.
function showChartGeneDetail(event) {
  byId("gene-detail").textContent = event.currentTarget.dataset.detail;
}
// Re-query the same selection with a smaller flanking window.
function zoomInRegion() {
  byId("window").selectedIndex = Math.max(0, byId("window").selectedIndex - 1);
  applyFilters();
}
// Re-query the same selection with a larger flanking window.
function zoomOutRegion() {
  byId("window").selectedIndex = Math.min(
    byId("window").options.length - 1,
    byId("window").selectedIndex + 1,
  );
  applyFilters();
}
// Apply basic filters to the selected gene or coordinate interval.
function applyFilters(event) {
  if (event) event.preventDefault();
  if (selectedGene) loadRegion(selectedGene, 1);
}
// Restore unfiltered association results without discarding the genomic selection.
function resetFilters() {
  byId("trait").value = "";
  byId("maximum-p").value = "1";
  byId("position").value = "";
  byId("include-repeats").checked = false;
  sortOrder = "p_asc";
  applyFilters();
}
// Toggle server-side P-value ordering across the complete filtered result set.
function toggleSort() {
  sortOrder = sortOrder === "p_asc" ? "p_desc" : "p_asc";
  applyFilters();
}
// Open the preceding association results page.
function previousPage() {
  if (selectedGene) loadRegion(selectedGene, currentPage - 1);
}
// Open the next association results page.
function nextPage() {
  if (selectedGene) loadRegion(selectedGene, currentPage + 1);
}
// Navigate to a requested page after validating its whole-number index.
function goToPage(event) {
  event.preventDefault();
  const page = Number(byId("page-number").value);
  if (Number.isInteger(page) && page > 0 && selectedGene)
    loadRegion(selectedGene, page);
}
// Populate a documented example and perform the same user-facing search.
function chooseExample(event) {
  byId("query").value = event.currentTarget.dataset.query;
  searchGenes();
}
// Restore bookmarked filters and selection when loading a page or using browser history.
function restoreBookmark() {
  const query = new URLSearchParams(location.hash.slice(1));
  ++latestRequest;
  selectedGene = query.get("gene") || "";
  byId("query").value = selectedGene;
  const windowValue = query.get("window") || "25000";
  byId("window").value = Array.from(byId("window").options).some(matchesWindow)
    ? windowValue
    : "25000";
  // Preserve only one of the supported sidebar window choices.
  function matchesWindow(option) {
    return option.value === windowValue;
  }
  byId("trait").value = query.get("trait") || "";
  byId("maximum-p").value = query.get("maximum_p") || "1";
  byId("position").value = query.get("position") || "";
  byId("include-repeats").checked = query.get("include_repeats") === "true";
  sortOrder = query.get("sort") === "p_desc" ? "p_desc" : "p_asc";
  currentPage = Math.max(1, Number(query.get("page")) || 1);
  showToolPage();
  if (selectedGene) loadRegion(selectedGene, currentPage);
  else {
    regionData = null;
    byId("selection").hidden = true;
    for (const panel of document.querySelectorAll(".requires-selection"))
      panel.hidden = true;
    for (const panel of document.querySelectorAll(".empty-selection"))
      panel.hidden = false;
    updateToolLinks();
  }
}
// Populate homepage statistics from the serving databases rather than static marketing numbers.
async function loadOverview() {
  try {
    const data = await fetchJson("/api/overview");
    byId("study-count").textContent = data.study_contexts.toLocaleString();
    byId("publication-count").textContent =
      `Across ${data.publications.toLocaleString()} publications. One publication may contain several study contexts.`;
    byId("trait-count").textContent = data.traits.toLocaleString();
    byId("cell-count").textContent = data.atlas_cells.toLocaleString();
    byId("cell-count-note").textContent = data.cell_count_note;
    byId("atlas-datasets").replaceChildren();
    for (const dataset of data.atlas_datasets) {
      const item = document.createElement("li");
      item.textContent = `${dataset.citation}: ${dataset.n_cells.toLocaleString()} · ${dataset.organ}`;
      byId("atlas-datasets").append(item);
    }
  } catch (error) {
    byId("overview-message").textContent =
      "Database counts could not be loaded. Please reload the page.";
  }
}
// Explain the planned homolog search without inventing cross-species results.
function submitHomolog(event) {
  event.preventDefault();
  byId("homolog-message").textContent =
    "Homolog searches will be available after a validated cross-species mapping source is integrated.";
}
// Dismiss suggestions when the user clicks outside the shared search box.
function dismissSuggestions(event) {
  if (!event.target.closest(".suggestion-anchor")) closeSuggestions();
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
  const width = Math.max(560, svg.parentElement.clientWidth);
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
    const fraction = (position - data.start) / Math.max(1, data.end - data.start);
    return leftMargin + fraction * (width - leftMargin - rightMargin);
  }

  // Convert an association strength into a vertical chart coordinate.
  function y(value) {
    return topMargin + plotHeight - (value / maximum) * plotHeight;
  }

  const bandLeft = x(data.gene.start);
  const bandRight = x(data.gene.end);
  let markup = `<rect x="${bandLeft}" y="${topMargin}" width="${Math.max(2, bandRight - bandLeft)}" height="${plotHeight}" fill="#0072b2" fill-opacity=".09"/><line x1="${bandLeft}" x2="${bandLeft}" y1="${topMargin}" y2="225" stroke="#0072b2" stroke-dasharray="4 4"/><line x1="${bandRight}" x2="${bandRight}" y1="${topMargin}" y2="225" stroke="#0072b2" stroke-dasharray="4 4"/>`;
  for (let index = 0; index <= 4; index++) {
    const value = (maximum * index) / 4;
    markup += `
      <line x1="${leftMargin}" x2="${width - rightMargin}" y1="${y(value)}" y2="${y(value)}" stroke="#e1e7ef"/>
      <text x="${leftMargin - 8}" y="${y(value) + 4}" text-anchor="end" font-size="14" fill="#647086">${value.toFixed(1)}</text>
    `;
  }
  markup +=
    '<text transform="translate(16,135) rotate(-90)" font-size="14" fill="#647086">−log₁₀(reported p)</text>';

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
  const tickCount = width < 800 ? 2 : 4;
  for (let index = 0; index <= tickCount; index++) {
    const position = Math.round(
      data.start + ((data.end - data.start) * index) / tickCount,
    );
    markup += `<text x="${x(position)}" y="210" text-anchor="${index === tickCount ? "end" : "middle"}" font-size="14" fill="#647086">${position.toLocaleString()}</text>`;
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
        ${compact ? "" : `<text x="${Math.min(left, width - rightMargin - 150)}" y="${centerY + 21}" font-size="14" fill="${color}">${escapeHtml(gene.symbol || gene.id)}</text>`}
        <title>${escapeHtml(label)}</title></g>`;
  }
  if (compact) {
    const selected = data.gene;
    const labelX = Math.max(
      leftMargin + 85,
      Math.min(width - rightMargin - 85, x((selected.start + selected.end) / 2)),
    );
    markup += `<line x1="${x((selected.start + selected.end) / 2)}" x2="${labelX}" y1="${250 + 3 * 22 + 12}" y2="345" stroke="#005a8d"/><text x="${labelX}" y="363" text-anchor="middle" font-size="14" fill="#005a8d">${escapeHtml(selected.symbol || selected.id)} · selected gene</text>`;
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

for (const link of document.querySelectorAll("[data-route]"))
  link.onclick = navigateTool;
for (const example of document.querySelectorAll(".example"))
  example.onclick = chooseExample;
byId("search-form").onsubmit = searchGenes;
byId("query").oninput = queueSuggestions;
byId("query").onkeydown = navigateSuggestions;
byId("window").onchange = applyFilters;
byId("filter-form").onsubmit = applyFilters;
byId("reset-filters").onclick = resetFilters;
byId("include-repeats").onchange = applyFilters;
byId("sort-p").onclick = toggleSort;
byId("prev").onclick = previousPage;
byId("next").onclick = nextPage;
byId("page-form").onsubmit = goToPage;
byId("zoom-in").onclick = zoomInRegion;
byId("zoom-out").onclick = zoomOutRegion;
byId("neighbor-search").oninput = filterNeighborTable;
byId("homolog-form").onsubmit = submitHomolog;
document.addEventListener("click", dismissSuggestions);
window.addEventListener("popstate", restoreBookmark);
window.addEventListener("hashchange", restoreBookmark);
// Keep the chart labels readable when the page width changes.
function resizeChart() {
  if (currentTool() === "gwas" && regionData) drawNeighborhood(regionData);
}
window.addEventListener("resize", resizeChart);
let atlasData = null;
let atlasRequest = 0;
let atlasGenes = [];
// Load only the validated full Wang atlas and retain the shared selected gene.
async function loadWangAtlas() {
  const request = ++atlasRequest;
  byId("atlas-content").hidden = true;
  if (!regionData || !regionData.gene.id) {
    byId("atlas-message").textContent =
      "Select a gene in the sidebar or Locus and Genes to explore its expression.";
    return;
  }
  if (!atlasGenes.includes(selectedGene)) atlasGenes = [selectedGene];
  byId("atlas-message").textContent = "Loading the eight-organ Wang atlas…";
  try {
    const data = await fetchJson(
      "/api/expression?" +
        new URLSearchParams({
          dataset: "wang_2025_eight_organ_atlas",
          genes: atlasGenes.join(","),
        }),
    );
    if (request !== atlasRequest) return;
    if (!data.available) {
      byId("atlas-message").textContent =
        "The full Wang atlas is not installed and validated in this environment yet.";
      return;
    }
    atlasData = data;
    const previous = byId("atlas-organ").value;
    byId("atlas-organ").replaceChildren(new Option("All eight organs", ""));
    for (const organ of data.dataset.organs)
      byId("atlas-organ").add(new Option(organ, organ));
    byId("atlas-organ").value = data.dataset.organs.includes(previous)
      ? previous
      : "Crown root";
    byId("atlas-compare").replaceChildren(new Option("Choose a neighbouring gene", ""));
    for (const gene of neighboringGenes)
      if (!atlasGenes.includes(gene.id))
        byId("atlas-compare").add(
          new Option(`${gene.symbol || gene.id} · ${gene.id}`, gene.id),
        );
    byId("atlas-source").innerHTML =
      `${escapeHtml(data.dataset.citation)} · ${data.dataset.n_cells.toLocaleString()} deposited cells · 8 organs · 16 samples. <a href="${escapeHtml(data.dataset.publication_url)}" target="_blank" rel="noopener noreferrer">Publication</a> · <a href="${escapeHtml(data.dataset.repository_url)}" target="_blank" rel="noopener noreferrer">Original GEO data</a>`;
    byId("atlas-export").href =
      "/api/expression/export?" +
      new URLSearchParams({ dataset: data.dataset.id, genes: atlasGenes.join(",") });
    byId("atlas-message").textContent = "";
    byId("atlas-content").hidden = false;
    renderAtlas();
  } catch (error) {
    if (request === atlasRequest) byId("atlas-message").textContent = error.message;
  }
}
// Select one organ without pooling measurements from other studies or organs.
function renderAtlas() {
  if (!atlasData) return;
  const organ = byId("atlas-organ").value;
  const contexts = atlasData.cell_groups.filter(matchesAtlasOrgan);
  // Retain all original organ-specific contexts selected by the organ filter.
  function matchesAtlasOrgan(group) {
    return !organ || group.organ === organ;
  }
  byId("atlas-genes").innerHTML = atlasGenes.map(renderAtlasGeneChip).join("");
  for (const button of byId("atlas-genes").querySelectorAll("[data-remove-gene]"))
    button.onclick = removeAtlasGene;
  let maximum = 0;
  const lookup = new Map();
  const values = [];
  for (const gene of atlasData.genes) {
    lookup.set(gene.gene_id, new Map());
    if (gene.status !== "measured") continue;
    for (const group of gene.groups) {
      lookup.get(gene.gene_id).set(group.cell_type, group);
      if (!organ || group.organ === organ) {
        maximum = Math.max(maximum, group.mean_log);
        values.push({ ...group, gene_id: gene.gene_id });
      }
    }
  }
  const left = organ ? 260 : 390;
  const width = Math.max(
    560,
    byId("atlas-chart").parentElement.clientWidth,
    left + atlasGenes.length * 125 + 30,
  );
  byId("atlas-chart").style.width = width + "px";
  const top = 65,
    spacing = (width - left - 30) / Math.max(1, atlasGenes.length);
  let markup = `<defs><linearGradient id="atlas-scale"><stop stop-color="#edf5fa"/><stop offset="1" stop-color="#0072b2"/></linearGradient></defs>`;
  for (let index = 0; index < atlasGenes.length; index++)
    markup += `<text x="${left + spacing * (index + 0.5)}" y="25" text-anchor="middle" font-size="14" fill="#233247">${escapeHtml(atlasGenes[index])}</text>`;
  for (let index = 0; index < contexts.length; index++) {
    const context = contexts[index];
    const y = top + index * 30;
    if (index % 2 === 0)
      markup += `<rect x="0" y="${y - 14}" width="${width}" height="30" fill="#f3f7f9"/>`;
    const label = organ
      ? context.source_cell_type
      : `${context.organ} · ${context.source_cell_type}`;
    markup += `<text x="${left - 12}" y="${y + 4}" text-anchor="end" font-size="14" fill="#233247">${escapeHtml(label)}</text>`;
    for (let column = 0; column < atlasGenes.length; column++) {
      const gene = atlasGenes[column],
        group = lookup.get(gene)?.get(context.cell_type),
        x = left + spacing * (column + 0.5);
      if (!group) {
        markup += `<text x="${x}" y="${y + 4}" text-anchor="middle" font-size="14" fill="#586779">N/A</text>`;
        continue;
      }
      const ratio = maximum ? group.mean_log / maximum : 0;
      const red = Math.round(237 * (1 - ratio)),
        green = Math.round(245 * (1 - ratio) + 114 * ratio),
        blue = Math.round(250 * (1 - ratio) + 178 * ratio);
      const radius = group.percent_expressing
        ? 3 + 9 * Math.sqrt(group.percent_expressing / 100)
        : 1.5;
      const text = `${gene} · ${context.organ} · ${context.source_cell_type}: mean log1p CP10K ${group.mean_log.toFixed(5)}, ${group.percent_expressing.toFixed(2)}% expressing (${group.n_expressing}/${group.n_cells} cells)`;
      markup += `<circle cx="${x}" cy="${y}" r="${radius}" fill="rgb(${red},${green},${blue})" stroke="#62798b" stroke-width=".7" role="img" aria-label="${escapeHtml(text)}"><title>${escapeHtml(text)}</title></circle>`;
    }
  }
  const bottom = top + contexts.length * 30 + 28;
  markup += `<rect x="${left}" y="${bottom}" width="150" height="12" fill="url(#atlas-scale)"/><text x="${left}" y="${bottom + 30}" font-size="14" fill="#586779">0</text><text x="${left + 150}" y="${bottom + 30}" text-anchor="end" font-size="14" fill="#586779">${maximum.toFixed(3)}</text><text x="${left}" y="${bottom + 52}" font-size="14" fill="#586779">Mean log1p(CP10K)</text>`;
  for (const [index, percent] of [0, 25, 50, 100].entries()) {
    const x = 45 + index * 55;
    const radius = percent ? 3 + 9 * Math.sqrt(percent / 100) : 1.5;
    markup += `<circle cx="${x}" cy="${bottom + 6}" r="${radius}" fill="#0072b2" stroke="#62798b" stroke-width=".7"/><text x="${x}" y="${bottom + 30}" text-anchor="middle" font-size="14" fill="#586779">${percent}%</text>`;
  }
  markup += `<text x="25" y="${bottom + 52}" font-size="14" fill="#586779">Expressing cells</text>`;
  byId("atlas-chart").setAttribute("viewBox", `0 0 ${width} ${bottom + 70}`);
  byId("atlas-chart").innerHTML = markup;
  byId("atlas-values").innerHTML = values.length
    ? values.map(renderAtlasValue).join("")
    : '<tr><td colspan="7">No mapped measurements are available for the selected genes. An absent feature is not zero expression.</td></tr>';
  const selected = atlasData.genes.find(isSelectedAtlasGene);
  // Identify the selected sidebar gene independently of comparison-gene ordering.
  function isSelectedAtlasGene(gene) {
    return gene.gene_id === selectedGene;
  }
  const groups = (selected?.groups || []).filter(matchesAtlasOrgan);
  byId("atlas-summary").textContent = groups.length
    ? `Highest mean expression ${organ ? "in " + organ : "across the deposited atlas"}: ${groups[0].source_cell_type} (${groups[0].organ}), ${groups[0].mean_log.toFixed(3)} mean log1p(CP10K), ${groups[0].percent_expressing.toFixed(1)}% expressing. Wang et al. (2025).`
    : `${selectedGene} is not mapped to an available feature in this deposited count assay; it is not treated as zero expression.`;
}
// Present each expression comparison as an explicitly removable gene choice.
function renderAtlasGeneChip(gene) {
  return `<span class="atlas-chip">${escapeHtml(gene)}${gene === selectedGene ? " · selected" : ` <button type="button" class="text-button" data-remove-gene="${escapeHtml(gene)}" aria-label="Remove ${escapeHtml(gene)} from comparison">Remove</button>`}</span>`;
}
// Retain exact numeric expression evidence as an accessible alternative to the dot plot.
function renderAtlasValue(row) {
  return `<tr><td>${escapeHtml(row.gene_id)}</td><td>${escapeHtml(row.organ)}</td><td>${escapeHtml(row.source_cell_type)}</td><td>${row.mean_log.toFixed(6)}</td><td>${row.percent_expressing.toFixed(2)}%</td><td>${row.n_expressing.toLocaleString()} / ${row.n_cells.toLocaleString()}</td><td>${row.rank}</td></tr>`;
}
// Add a neighbouring gene for comparison within the same validated study.
function addAtlasGene() {
  const gene = byId("atlas-compare").value;
  if (!gene || atlasGenes.includes(gene)) return;
  if (atlasGenes.length >= 6) {
    byId("atlas-message").textContent =
      "Compare up to six genes at a time to keep the plot readable.";
    return;
  }
  atlasGenes.push(gene);
  loadWangAtlas();
}
// Remove only a comparison gene while preserving the shared selected gene.
function removeAtlasGene(event) {
  // Preserve every comparison except the gene explicitly removed by the user.
  function keepAtlasGene(gene) {
    return gene !== event.currentTarget.dataset.removeGene;
  }
  atlasGenes = atlasGenes.filter(keepAtlasGene);
  loadWangAtlas();
}
byId("atlas-organ").onchange = renderAtlas;
byId("atlas-add").onclick = addAtlasGene;

restoreBookmark();
loadOverview();
