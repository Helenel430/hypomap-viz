// Port of ../../gene-explorer/app.R (R/Shiny) to plain JS.
// Loads the precomputed data/gene_stats.json (see scripts/export_gene_stats.R)
// and reproduces the gene lookup, suggestions, summary, and per-cell-type bars.

const app = document.getElementById("app");

async function main() {
  const stats = await fetch("/data/gene_stats.json").then((r) => r.json());
  render(stats);
}

function render(stats) {
  const genes = stats.genes;
  const genesLc = genes.map((g) => g.toLowerCase());
  const geneIndexByLc = new Map(genesLc.map((g, i) => [g, i]));
  const cts = stats.celltypes;
  const total = stats.totalCells;
  const nExpressedGlobal = stats.global.nExpressing.filter((n) => n > 0).length;

  // 99th percentile of mean_expressing (expressing cells) across gene x celltype,
  // used as the fixed scale for the "Mean" bar charts — mirrors app.R's mean_scale.
  const meanScale = percentile99(
    stats.byct.meanExpressing.flatMap((row, i) =>
      row.filter((_, ct) => stats.byct.nExpressing[i][ct] > 0)
    )
  );

  app.innerHTML = `
    <div class="sidebar">
      <h1>Romanov_10x — Gene Expression Lookup</h1>
      <label for="gene_in">Explore a gene:</label>
      <input type="text" id="gene_in" placeholder="Start typing a gene (case-insensitive)..." value="Pomc" autocomplete="off" />
      <div class="suggestions" id="suggestions"></div>
      <button class="primary" id="go">Look up</button>
      <p class="help">Search for a gene (e.g. igf) and click Look up or press Enter.</p>
      <p class="help">
        <b>Data:</b> the HypoMap P23-only subset of <b>Romanov_10x</b>, a log-normalized
        single-cell RNA dataset of developing mouse hypothalamus
        (${total.toLocaleString()} cells, ${cts.length} cell types).
      </p>
    </div>
    <div class="main">
      <h2 id="gtitle"></h2>
      <div id="globalbox"></div>
      <hr />
      <h3>By cell type</h3>
      <div id="cttable"></div>
      <h3>% expressing</h3>
      <div id="pctbars"></div>
      <h3>Mean (expressing cells)</h3>
      <div id="meanbars"></div>
    </div>
  `;

  const geneInput = document.getElementById("gene_in");
  const suggestionsEl = document.getElementById("suggestions");
  const goBtn = document.getElementById("go");

  function matchGenes(q, n = 12) {
    q = q.trim().toLowerCase();
    if (!q) return [];
    const exact = genes.filter((_, i) => genesLc[i] === q);
    const pre = genes.filter((_, i) => genesLc[i] !== q && genesLc[i].startsWith(q));
    const sub = genes.filter((_, i) => !genesLc[i].startsWith(q) && genesLc[i].includes(q));
    pre.sort();
    sub.sort();
    return Array.from(new Set([...exact, ...pre, ...sub])).slice(0, n);
  }

  function suggestText(q) {
    const m = matchGenes(q, 6);
    if (m.length === 0) return "No similar symbols found.";
    return `Did you mean: ${m.join(", ")}?`;
  }

  function renderSuggestions() {
    const q = geneInput.value.trim();
    if (q.length < 2) {
      suggestionsEl.innerHTML = "";
      return;
    }
    if (geneIndexByLc.has(q.toLowerCase())) {
      suggestionsEl.innerHTML = "";
      return;
    }
    const m = matchGenes(q, 12);
    if (m.length === 0) {
      suggestionsEl.innerHTML = `<div class="empty">no matches</div>`;
      return;
    }
    suggestionsEl.innerHTML = m
      .map((g) => `<a data-gene="${g}">${g}</a>`)
      .join("");
    suggestionsEl.querySelectorAll("a").forEach((a) => {
      a.addEventListener("click", (e) => {
        e.preventDefault();
        geneInput.value = a.dataset.gene;
        suggestionsEl.innerHTML = "";
        lookup(a.dataset.gene);
      });
    });
  }

  function lookup(query) {
    const q = query.trim();
    const idx = geneIndexByLc.get(q.toLowerCase());
    if (idx === undefined) {
      showNotFound(q);
      return;
    }
    showGene(idx);
  }

  function showNotFound(q) {
    document.getElementById("gtitle").textContent = "";
    document.getElementById("globalbox").innerHTML =
      `<p class="placeholder">Gene '${escapeHtml(q)}' not found — try one of the suggestions. ${escapeHtml(suggestText(q))}</p>`;
    document.getElementById("cttable").innerHTML = "";
    document.getElementById("pctbars").innerHTML = "";
    document.getElementById("meanbars").innerHTML = "";
  }

  function showGene(idx) {
    const gname = genes[idx];
    document.getElementById("gtitle").textContent = gname;

    const nExp = stats.global.nExpressing[idx];
    const pctExp = stats.global.pctExpressing[idx];
    const meanExp = stats.global.meanExpressing[idx];
    const pctile = stats.global.percentile[idx];
    const ncells = total.toLocaleString();

    if (nExp === 0) {
      document.getElementById("globalbox").innerHTML =
        `<ul class="summary"><li><i>${gname}</i> is not detected in any of the ${ncells} P23 cells in this dataset.</li></ul>`;
    } else {
      document.getElementById("globalbox").innerHTML = `
        <ul class="summary">
          <li><i>${gname}</i> is expressed in <b>${pctExp.toFixed(2)}%</b> of cells — it is detected in
              ${nExp.toLocaleString()} of the ${ncells} P23 cells in this dataset.</li>
          <li>Averaged over only the cells that express it, its mean (log-normalized)
              expression is <b>${meanExp.toFixed(3)}</b>.</li>
          <li>Compared with the ${nExpressedGlobal.toLocaleString()} genes that are detectably expressed
              in this cohort, and ranked by that same mean-when-expressed level, <i>${gname}</i> sits in the
              <b>${pctile.toFixed(1)}th percentile</b> — its expression, when on, is higher than about
              ${pctile.toFixed(0)}% of those genes.</li>
        </ul>`;
    }

    const rows = cts.map((ct, ctIdx) => ({
      cellType: ct,
      nCells: stats.cellCounts[ctIdx],
      nExpressing: stats.byct.nExpressing[idx][ctIdx],
      pctExpressing: stats.byct.pctExpressing[idx][ctIdx],
      meanExpressing: stats.byct.meanExpressing[idx][ctIdx],
      percentile: stats.byct.percentile[idx][ctIdx],
    }));

    document.getElementById("cttable").innerHTML = renderTable(rows);
    document.getElementById("pctbars").innerHTML = renderBars(
      rows.map((r) => ({ label: r.cellType, value: r.pctExpressing, count: r.nExpressing })),
      100,
      (v) => `${v.toFixed(1)}%`,
      "#3aa37a"
    );
    document.getElementById("meanbars").innerHTML = renderBars(
      rows.map((r) => ({ label: r.cellType, value: r.meanExpressing, count: r.nExpressing })),
      meanScale,
      (v) => v.toFixed(2),
      "#3a6ea3"
    );
  }

  function renderTable(rows) {
    const sorted = [...rows].sort((a, b) => b.pctExpressing - a.pctExpressing);
    const body = sorted
      .map((r) => {
        const top5 = r.percentile >= 95 ? "Yes" : "";
        return `<tr>
          <td>${r.cellType}</td>
          <td>${r.nCells.toLocaleString()}</td>
          <td>${r.nExpressing.toLocaleString()}</td>
          <td>${r.pctExpressing.toFixed(2)}</td>
          <td>${r.meanExpressing.toFixed(2)}</td>
          <td>${r.percentile != null && isFinite(r.percentile) ? r.percentile.toFixed(2) : ""}</td>
          <td>${top5}</td>
        </tr>`;
      })
      .join("");
    return `<table class="ct-table">
      <thead><tr>
        <th>Cell type</th><th>#Cells</th><th>#Expressing</th><th>%Expressing</th>
        <th>Mean expr. (expressing cells)</th><th>Percentile (within cell type)</th><th>Top 5% (in cell type)</th>
      </tr></thead>
      <tbody>${body}</tbody>
    </table>`;
  }

  function renderBars(items, maxval, fmt, accent) {
    const filtered = items.filter((it) => it.value != null && !isNaN(it.value));
    if (filtered.length === 0 || Math.max(...filtered.map((it) => it.value)) <= 0) {
      return `<div class="placeholder">Not expressed in any cell type.</div>`;
    }
    filtered.sort((a, b) => b.value - a.value);
    return filtered
      .map((it) => {
        const over = it.value > maxval;
        const w = Math.min(100, (100 * it.value) / maxval);
        return `<div class="bar-row">
          <div class="bar-label">${it.label}</div>
          <div class="bar-line">
            <div class="bar-track">
              <div class="bar-fill" style="width:${w.toFixed(1)}%;background:${accent}"></div>
            </div>
            <div class="bar-value" style="color:${over ? accent : "#333"}">${over ? "&#9656; " : ""}${fmt(it.value)}</div>
            <div class="bar-count">${it.count.toLocaleString()}</div>
          </div>
        </div>`;
      })
      .join("");
  }

  geneInput.addEventListener("input", renderSuggestions);
  geneInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      suggestionsEl.innerHTML = "";
      lookup(geneInput.value);
    }
  });
  goBtn.addEventListener("click", () => {
    suggestionsEl.innerHTML = "";
    lookup(geneInput.value);
  });

  lookup("Pomc");
}

function percentile99(values) {
  const sorted = [...values].sort((a, b) => a - b);
  const idx = Math.floor(0.99 * (sorted.length - 1));
  return sorted[idx];
}

function escapeHtml(s) {
  const div = document.createElement("div");
  div.textContent = s;
  return div.innerHTML;
}

main();
