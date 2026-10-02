/*
 * perfcharts.js - small interactive SVG charts for the "Performance" pages.
 *
 * Self-contained (no dependencies, no network besides the chart's own JSON
 * file next to the page). A chart is a Markdown link with the attr_list class
 * `perf-chart`, alone in its paragraph:
 *
 *   [Figure 1 data (JSON)](data/<name>.json){ .perf-chart }
 *
 * whose JSON is written by baysor-benchmarks/docs_figures/make_figures.py.
 * The script replaces the paragraph with the chart: panels of log/linear/band
 * axes, line and marker series, reference lines, a legend, a hover/keyboard
 * tooltip and a data table. Without JavaScript the link to the data remains.
 * Colours come from CSS variables (perfcharts.css), so the charts follow the
 * Material colour scheme and its toggle without redrawing.
 *
 * It also makes a table sortable by column when the paragraph right after it
 * carries the class `perf-sortable` (`{: .perf-sortable }`), as do the charts'
 * data tables; numbers with units (s, min, h, B, KiB, MiB, GiB, %) and
 * thousands separators sort by value.
 */
(function () {
  "use strict";
  var NS = "http://www.w3.org/2000/svg";
  var DASH = { solid: null, dashed: "6 4", dotted: "1.5 3.5" };

  function el(tag, attrs, parent) {
    var e = document.createElementNS(NS, tag);
    for (var k in attrs) if (attrs[k] !== null && attrs[k] !== undefined) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function html(tag, cls, parent, text) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined) e.textContent = text;
    if (parent) parent.appendChild(e);
    return e;
  }
  function color(c) {
    return typeof c === "number" ? "var(--pc-c" + c + ")" : "var(--pc-" + (c || "ink2") + ")";
  }

  // ------------------------------------------------------------ scales ----
  function niceLinear(lo, hi, n) {
    var span = hi - lo || Math.abs(hi) || 1;
    var step = Math.pow(10, Math.floor(Math.log10(span / n)));
    var err = span / n / step;
    step *= err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1;
    var ticks = [];
    for (var v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) ticks.push(+v.toPrecision(12));
    return ticks;
  }
  function logTicks(lo, hi, px) {
    var all = [], ones = [];
    for (var e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) {
      [1, 2, 5].forEach(function (d) {
        var v = d * Math.pow(10, e);
        if (v >= lo * 0.999 && v <= hi * 1.001) { all.push(v); if (d === 1) ones.push(v); }
      });
    }
    var maxTicks = Math.max(3, Math.floor(px / 45));
    if (all.length <= maxTicks) return all;
    var two = all.filter(function (v) { return !/^2/.test(v.toPrecision(1)); });
    return two.length <= maxTicks || ones.length < 2 ? two : ones;
  }
  function scale(axis, data, px0, px1) {
    if (axis.type === "band") {
      var cats = axis.categories, step = (px1 - px0) / cats.length;
      var f = function (v) { return px0 + step * (cats.indexOf(v) + 0.5); };
      f.ticks = cats; f.band = step;
      return f;
    }
    var lo = axis.min, hi = axis.max;
    var vals = data.filter(function (v) { return v !== null && isFinite(v) && (!axis.log || v > 0); });
    if (lo === undefined) lo = Math.min.apply(null, vals);
    if (hi === undefined) hi = Math.max.apply(null, vals);
    if (axis.log) {
      var pad = (Math.log10(hi) - Math.log10(lo) || 1) * 0.06;
      if (axis.min === undefined) lo = Math.pow(10, Math.log10(lo) - pad);
      if (axis.max === undefined) hi = Math.pow(10, Math.log10(hi) + pad);
      var a = Math.log10(lo), b = Math.log10(hi);
      var g = function (v) { return px0 + (Math.log10(v) - a) / (b - a) * (px1 - px0); };
      g.ticks = axis.ticks || logTicks(lo, hi, Math.abs(px1 - px0));
      return g;
    }
    if (axis.min === undefined) lo = Math.min(lo, axis.zero === false ? lo : 0);
    var p = (hi - lo) * 0.05;
    if (axis.max === undefined) hi += p;
    var t = axis.ticks || niceLinear(lo, hi, Math.max(3, Math.floor(Math.abs(px1 - px0) / 60)));
    if (axis.max === undefined) hi = Math.max(hi, t[t.length - 1]);
    var h = function (v) { return px0 + (v - lo) / (hi - lo) * (px1 - px0); };
    h.ticks = t;
    return h;
  }
  function fmtTick(v, fmt) {
    if (fmt === "plain") return v.toLocaleString("en-US", { maximumFractionDigits: 3 });
    if (fmt === "pct") return v + " %";
    var a = Math.abs(v);
    if (a >= 1e6) return +(v / 1e6).toPrecision(3) + "M";
    if (a >= 1e3) return +(v / 1e3).toPrecision(3) + "k";
    return String(+v.toPrecision(3));
  }

  // ----------------------------------------------------------- markers ----
  function marker(g, shape, x, y, r, s) {
    var m;
    if (shape === "square") {
      m = el("rect", { x: x - r * 0.9, y: y - r * 0.9, width: r * 1.8, height: r * 1.8 }, g);
    } else if (shape === "triangle") {
      m = el("path", { d: "M" + x + "," + (y - r * 1.2) + "L" + (x + r * 1.1) + "," + (y + r * 0.75) +
                       "L" + (x - r * 1.1) + "," + (y + r * 0.75) + "Z" }, g);
    } else if (shape === "diamond") {
      m = el("path", { d: "M" + x + "," + (y - r * 1.25) + "L" + (x + r * 1.25) + "," + y + "L" + x + "," +
                       (y + r * 1.25) + "L" + (x - r * 1.25) + "," + y + "Z" }, g);
    } else {
      m = el("circle", { cx: x, cy: y, r: r }, g);
    }
    var c = color(s.color);
    if (s.open) m.setAttribute("style", "fill:var(--pc-surface);stroke:" + c + ";stroke-width:1.6");
    else m.setAttribute("style", "fill:" + c + ";stroke:var(--pc-surface);stroke-width:1.5");
    return m;
  }
  function legendKey(item) {
    var svg = el("svg", { width: 26, height: 14, "aria-hidden": "true" });
    if (item.line) {
      el("line", { x1: 1, x2: 25, y1: 7, y2: 7, "stroke-width": 2, "stroke-dasharray": DASH[item.line],
                   style: "stroke:" + color(item.color) }, svg);
    }
    if (item.marker) marker(svg, item.marker, 13, 7, 4, item);
    if (item.rect) el("rect", { x: 4, y: 2, width: 18, height: 10, rx: 2, style: "fill:" + color(item.color) }, svg);
    return svg;
  }

  // ------------------------------------------------------------- panel ----
  function drawPanel(box, panel, spec, width, tip) {
    var wrap = html("div", "perf-panel", box);
    if (panel.title) html("div", "perf-panel-title", wrap, panel.title);
    var H = panel.height || spec.height || 300;
    var m = { l: panel.marginLeft || 58, r: panel.marginRight || 18, t: 10, b: 44 };
    var svg = el("svg", { width: width, height: H, role: "img", tabindex: 0, class: "perf-svg",
                          "aria-label": (panel.title || spec.title || "chart") +
                          ". Use the arrow keys to read the points, or open the data table below." }, wrap);
    var xs = [], ys = [];
    panel.series.forEach(function (s) { s.points.forEach(function (p) { xs.push(p.x); ys.push(p.y); }); });
    (panel.refs || []).forEach(function (r) {
      if (r.type === "line") r.points.forEach(function (p) { xs.push(p[0]); ys.push(p[1]); });
    });
    var X = scale(panel.x, xs, m.l, width - m.r), Y = scale(panel.y, ys, H - m.b, m.t);
    var gGrid = el("g", { class: "perf-grid" }, svg);
    // y grid + ticks
    Y.ticks.forEach(function (v) {
      var y = Y(v);
      if (panel.y.type !== "band") el("line", { x1: m.l, x2: width - m.r, y1: y, y2: y }, gGrid);
      var t = el("text", { x: m.l - 7, y: y, "text-anchor": "end", "dominant-baseline": "middle",
                           class: "perf-tick" }, svg);
      t.textContent = panel.y.type === "band" ? v : fmtTick(v, panel.y.fmt);
    });
    X.ticks.forEach(function (v) {
      var x = X(v);
      if (panel.x.type !== "band") el("line", { x1: x, x2: x, y1: m.t, y2: H - m.b }, gGrid);
      var t = el("text", { x: x, y: H - m.b + 16, "text-anchor": "middle", class: "perf-tick" }, svg);
      t.textContent = fmtTick(v, panel.x.fmt);
    });
    el("line", { x1: m.l, x2: width - m.r, y1: H - m.b, y2: H - m.b, class: "perf-axis" }, svg);
    if (panel.y.type !== "band") el("line", { x1: m.l, x2: m.l, y1: m.t, y2: H - m.b, class: "perf-axis" }, svg);
    var xl = el("text", { x: (m.l + width - m.r) / 2, y: H - 6, "text-anchor": "middle", class: "perf-label" }, svg);
    xl.textContent = panel.x.label || "";
    if (panel.y.label) {
      var yl = el("text", { transform: "translate(13," + (m.t + H - m.b) / 2 + ") rotate(-90)",
                            "text-anchor": "middle", class: "perf-label" }, svg);
      yl.textContent = panel.y.label;
    }
    // reference marks
    (panel.refs || []).forEach(function (r) {
      var g = el("g", { class: "perf-ref" }, svg), lx, ly;
      if (r.type === "vline") {
        lx = X(r.x);
        el("line", { x1: lx, x2: lx, y1: m.t, y2: H - m.b }, g);
        ly = m.t + 10;
      } else {
        var d = r.points.map(function (p, i) { return (i ? "L" : "M") + X(p[0]) + "," + Y(p[1]); }).join("");
        el("path", { d: d, "stroke-dasharray": DASH[r.dash || "dotted"] }, g);
        lx = X(r.labelAt ? r.labelAt[0] : r.points[1][0]);
        ly = Y(r.labelAt ? r.labelAt[1] : r.points[1][1]);
      }
      if (r.label) {
        var t = el("text", { x: lx + 4, y: ly, class: "perf-ref-label" }, g);
        t.textContent = r.label;
      }
    });
    (panel.links || []).forEach(function (k) {
      el("line", { x1: X(k.x1), x2: X(k.x2), y1: Y(k.y), y2: Y(k.y), class: "perf-link" }, svg);
    });
    // series
    var pts = [];
    panel.series.forEach(function (s) {
      var g = el("g", {}, svg);
      var P = s.points.filter(function (p) { return p.y !== null && p.x !== null; });
      if (s.line && P.length > 1) {
        el("path", { d: P.map(function (p, i) { return (i ? "L" : "M") + X(p.x) + "," + Y(p.y); }).join(""),
                     fill: "none", "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round",
                     "stroke-dasharray": DASH[s.line], style: "stroke:" + color(s.color) }, g);
      }
      P.forEach(function (p) {
        var px = X(p.x), py = Y(p.y);
        if (s.marker) marker(g, s.marker, px, py, s.size || 4, s);
        pts.push({ px: px, py: py, p: p, s: s });
      });
    });
    (panel.labels || []).forEach(function (l) {
      var x = X(l.x), dx = l.dx || 0;
      var t = el("text", { x: x + dx, y: Y(l.y) + (l.dy || 0), "text-anchor": "start", class: "perf-annot" }, svg);
      t.textContent = l.text;
      if (x + dx + t.getComputedTextLength() > width - m.r) {  // keep the label inside the plot
        t.setAttribute("x", x - Math.abs(dx));
        t.setAttribute("text-anchor", "end");
      }
    });
    // hover layer
    var hl = el("circle", { r: 7, class: "perf-hover", visibility: "hidden" }, svg);
    var cross = el("line", { y1: m.t, y2: H - m.b, class: "perf-cross", visibility: "hidden" }, svg);
    var order = pts.slice().sort(function (a, b) { return a.px - b.px || a.py - b.py; });
    var current = -1;
    function show(hit, clientX, clientY) {
      if (!hit) { hide(); return; }
      hl.setAttribute("cx", hit.px); hl.setAttribute("cy", hit.py);
      hl.setAttribute("visibility", "visible");
      var rows, title;
      if (panel.hover === "x") {
        cross.setAttribute("x1", hit.px); cross.setAttribute("x2", hit.px);
        cross.setAttribute("visibility", "visible");
        title = hit.p.xTitle;
        rows = pts.filter(function (q) { return q.p.x === hit.p.x; })
                  .sort(function (a, b) { return a.py - b.py; })
                  .map(function (q) { return { v: q.p.yText, l: q.p.rowLabel || q.s.name, s: q.s }; });
      } else {
        title = hit.p.title;
        rows = (hit.p.rows || []).map(function (r) { return { v: r[0], l: r[1] }; });
        rows[0] && (rows[0].s = hit.s);
      }
      tip.show(title, rows, svg, hit.px, hit.py);
    }
    function hide() {
      hl.setAttribute("visibility", "hidden");
      cross.setAttribute("visibility", "hidden");
      tip.hide();
    }
    function nearest(ev) {
      var r = svg.getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top, best = null, bd = Infinity;
      pts.forEach(function (q) {
        var d = panel.hover === "x" ? Math.abs(q.px - x) + Math.abs(q.py - y) * 0.05 : Math.hypot(q.px - x, q.py - y);
        if (d < bd) { bd = d; best = q; }
      });
      return bd <= (panel.hover === "x" ? 40 : 28) ? best : null;
    }
    svg.addEventListener("pointermove", function (ev) { var h = nearest(ev); current = order.indexOf(h); show(h); });
    svg.addEventListener("pointerleave", hide);
    svg.addEventListener("blur", hide);
    svg.addEventListener("keydown", function (ev) {
      if (ev.key === "ArrowRight" || ev.key === "ArrowDown") current = Math.min(order.length - 1, current + 1);
      else if (ev.key === "ArrowLeft" || ev.key === "ArrowUp") current = Math.max(0, current - 1);
      else if (ev.key === "Escape") { hide(); return; }
      else return;
      ev.preventDefault();
      show(order[current]);
    });
  }

  // ------------------------------------------------------------ tooltip ----
  function tooltip(root) {
    var t = html("div", "perf-tip", root);
    t.setAttribute("role", "status");
    return {
      show: function (title, rows, svg, px, py) {
        t.textContent = "";
        if (title) html("div", "perf-tip-title", t, title);
        rows.forEach(function (r) {
          var row = html("div", "perf-tip-row", t);
          if (r.s) {
            var k = legendKey({ line: r.s.line || "solid", color: r.s.color });
            k.setAttribute("class", "perf-tip-key");
            row.appendChild(k);
          }
          html("strong", null, row, r.v);
          html("span", null, row, " " + r.l);
        });
        t.style.display = "block";
        var rr = root.getBoundingClientRect(), sr = svg.getBoundingClientRect();
        var x = sr.left - rr.left + px + 14, y = sr.top - rr.top + py + 14;
        if (x + t.offsetWidth > rr.width) x = Math.max(0, sr.left - rr.left + px - t.offsetWidth - 14);
        t.style.left = x + "px";
        t.style.top = y + "px";
      },
      hide: function () { t.style.display = "none"; }
    };
  }

  // -------------------------------------------------------------- chart ----
  function dataTable(root, table) {
    var d = html("details", "perf-chart-table", root);
    html("summary", null, d, "Data table");
    var tb = html("table", null, html("div", "perf-table-wrap", d)), tr = html("tr", null, html("thead", null, tb));
    table.columns.forEach(function (c) { html("th", null, tr, c); });
    var body = html("tbody", null, tb);
    table.rows.forEach(function (r) {
      var row = html("tr", null, body);
      r.forEach(function (v) { html("td", null, row, v === null ? "—" : String(v)); });
    });
    makeSortable(tb);
  }

  function render(root, spec) {
    root.textContent = "";
    root.classList.add("perf-chart-ready");
    if (spec.legend && spec.legend.length) {
      var lg = html("div", "perf-legend", root);
      spec.legend.forEach(function (it) {
        var s = html("span", "perf-legend-item", lg);
        s.appendChild(legendKey(it));
        html("span", null, s, it.name);
      });
    }
    var box = html("div", "perf-panels", root);
    var tip = tooltip(root);
    var draw = function () {
      box.textContent = "";
      var W = box.clientWidth || 600, cols = spec.columns || 1;
      if (W < 620) cols = 1;
      box.style.gridTemplateColumns = "repeat(" + cols + ", minmax(0, 1fr))";
      var w = Math.floor((W - (cols - 1) * 16) / cols);
      spec.panels.forEach(function (p) { drawPanel(box, p, spec, w, tip); });
    };
    draw();
    if (spec.table) dataTable(root, spec.table);
    var lastW = box.clientWidth;
    if (window.ResizeObserver) {
      new ResizeObserver(function () {
        if (Math.abs(box.clientWidth - lastW) > 4) { lastW = box.clientWidth; draw(); }
      }).observe(box);
    }
  }

  function initCharts() {
    document.querySelectorAll("a.perf-chart[href$='.json']").forEach(function (a) {
      var root = document.createElement("div");
      root.className = "perf-chart";
      var p = a.parentNode;
      if (p.tagName === "P" && p.textContent.trim() === a.textContent.trim()) p.replaceWith(root);
      else p.insertBefore(root, a.nextSibling);
      root.appendChild(a);
      fetch(a.href).then(function (r) { return r.json(); }).then(function (spec) { render(root, spec); })
        .catch(function (e) { root.setAttribute("data-error", String(e)); });
    });
  }

  // ------------------------------------------------------- sortable tables ----
  var UNITS = { s: 1, min: 60, h: 3600, b: 1, kib: 1024, mib: 1048576, gib: 1073741824, "%": 1, "×": 1 };
  function cellValue(text) {
    var t = text.replace(/[,  ]/g, "").trim();
    var m = /^([-+−]?\d*\.?\d+(?:e[-+]?\d+)?)\s*([a-zA-Z%×]*)$/.exec(t);
    if (!m) return null;
    var u = m[2].toLowerCase();
    if (u && !(u in UNITS)) return null;
    return parseFloat(m[1].replace("−", "-")) * (u ? UNITS[u] : 1);
  }
  function makeSortable(table) {
    if (table.dataset.sortable) return;
    table.dataset.sortable = "1";
    var heads = table.querySelectorAll("thead th");
    var body = table.tBodies[0];
    if (!body) return;
    heads.forEach(function (th, ci) {
      var b = document.createElement("button");
      b.type = "button";
      b.className = "perf-sort";
      while (th.firstChild) b.appendChild(th.firstChild);
      th.appendChild(b);
      th.setAttribute("aria-sort", "none");
      b.addEventListener("click", function () {
        var dir = th.getAttribute("aria-sort") === "ascending" ? -1 : 1;
        heads.forEach(function (h) { h.setAttribute("aria-sort", "none"); });
        th.setAttribute("aria-sort", dir === 1 ? "ascending" : "descending");
        var rows = Array.prototype.slice.call(body.rows);
        var keys = rows.map(function (r) { var c = r.cells[ci]; return c ? c.textContent.trim() : ""; });
        var nums = keys.map(cellValue);
        var numeric = keys.every(function (k, i) { return nums[i] !== null || k === "" || k === "—"; }) &&
                      nums.some(function (n) { return n !== null; });
        var idx = rows.map(function (_, i) { return i; });
        idx.sort(function (a, b) {
          if (numeric) {
            var x = nums[a], y = nums[b];
            if (x === null || y === null) return x === y ? a - b : x === null ? 1 : -1;  // blanks last
            return (x - y) * dir || a - b;
          }
          return keys[a].localeCompare(keys[b], undefined, { numeric: true }) * dir || a - b;
        });
        idx.forEach(function (i) { body.appendChild(rows[i]); });
      });
    });
  }
  function initTables() {
    document.querySelectorAll("p.perf-sortable").forEach(function (p) {
      var prev = p.previousElementSibling;
      var t = prev && (prev.tagName === "TABLE" ? prev : prev.querySelector("table"));
      if (t) makeSortable(t);
    });
  }

  function init() { initCharts(); initTables(); }
  if (window.document$ && window.document$.subscribe) window.document$.subscribe(init);
  else if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
