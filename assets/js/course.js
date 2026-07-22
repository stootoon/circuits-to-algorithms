/* ------------------------------------------------------------------
   Circuits to Algorithms — reading tools.

   Everything lives in localStorage, per browser. No account, no server.
   - progress:  which units you have marked complete
   - marks:     persistent highlights, each either a note or a flagged error
   Flagged errors can be pushed to GitHub as prefilled issues.
   ------------------------------------------------------------------ */
(function () {
  "use strict";

  var REPO = "stootoon/circuits-to-algorithms";
  var K_PROGRESS = "c2a.progress";
  var K_MARKS = "c2a.marks";
  var K_NAV = "c2a.nav";
  var K_NAVHIDDEN = "c2a.navHidden";
  var CTX = 40; // characters of context stored either side of a highlight

  var PATH = location.pathname.replace(/index\.html$/, "");

  /* ------------------------------- storage ------------------------------ */

  function load(key, dflt) {
    try { return JSON.parse(localStorage.getItem(key)) || dflt; }
    catch (e) { return dflt; }
  }
  function save(key, val) {
    try { localStorage.setItem(key, JSON.stringify(val)); return true; }
    catch (e) {
      console.warn("[c2a] could not save — storage full or blocked", e);
      return false;
    }
  }
  var getProgress = function () { return load(K_PROGRESS, {}); };
  var getMarks = function () { return load(K_MARKS, {}); };
  function marksHere() { return getMarks()[PATH] || []; }
  function setMarksHere(list) {
    var all = getMarks();
    if (list.length) all[PATH] = list; else delete all[PATH];
    save(K_MARKS, all);
  }
  function uid() { return Date.now().toString(36) + Math.random().toString(36).slice(2, 7); }

  function pageTitle() {
    var h = document.querySelector(".main-content h1");
    return (h ? h.textContent : document.title).replace(/\s+/g, " ").trim();
  }

  /* --------------------------- text anchoring --------------------------- */
  /* A highlight is stored as a text-quote anchor (exact + surrounding
     context) rather than a DOM path, so it survives edits to the page that
     do not touch the quoted sentence. MathJax output is skipped so that the
     flattened text is identical when saving and when restoring.            */

  function flatten(root) {
    var walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, {
      acceptNode: function (n) {
        for (var p = n.parentNode; p && p !== root; p = p.parentNode) {
          var tag = p.nodeName;
          if (tag === "SCRIPT" || tag === "STYLE" || tag === "MJX-CONTAINER") {
            return NodeFilter.FILTER_REJECT;
          }
          if (p.classList && p.classList.contains("c2a-ui")) return NodeFilter.FILTER_REJECT;
        }
        return NodeFilter.FILTER_ACCEPT;
      }
    });
    var text = "", map = [], n;
    while ((n = walker.nextNode())) {
      map.push({ node: n, start: text.length, len: n.nodeValue.length });
      text += n.nodeValue;
    }
    return { text: text, map: map };
  }

  function offsetOf(flat, node, off) {
    for (var i = 0; i < flat.map.length; i++) {
      if (flat.map[i].node === node) return flat.map[i].start + off;
    }
    return -1;
  }

  function commonPrefix(a, b) {
    var i = 0; while (i < a.length && i < b.length && a[i] === b[i]) i++; return i;
  }
  function commonSuffix(a, b) {
    var i = 0; while (i < a.length && i < b.length && a[a.length - 1 - i] === b[b.length - 1 - i]) i++; return i;
  }

  function locate(flat, mark) {
    if (!mark.exact) return null;
    var best = -1, bestScore = -1, i = -1;
    while ((i = flat.text.indexOf(mark.exact, i + 1)) !== -1) {
      var pre = flat.text.slice(Math.max(0, i - CTX), i);
      var suf = flat.text.slice(i + mark.exact.length, i + mark.exact.length + CTX);
      var score = commonSuffix(pre, mark.prefix || "") + commonPrefix(suf, mark.suffix || "");
      if (score > bestScore) { bestScore = score; best = i; }
    }
    return best < 0 ? null : { start: best, end: best + mark.exact.length };
  }

  function paint(flat, start, end, mark) {
    var painted = false;
    for (var i = 0; i < flat.map.length; i++) {
      var e = flat.map[i], nStart = e.start, nEnd = e.start + e.len;
      if (nEnd <= start || nStart >= end) continue;
      var from = Math.max(start, nStart) - nStart;
      var to = Math.min(end, nEnd) - nStart;
      var node = e.node;
      if (!node.parentNode) continue;
      if (to < node.nodeValue.length) node.splitText(to);
      if (from > 0) node = node.splitText(from);
      var m = document.createElement("mark");
      m.className = "c2a-mark c2a-" + mark.type;
      m.dataset.id = mark.id;
      m.id = "c2a-" + mark.id;
      node.parentNode.replaceChild(m, node);
      m.appendChild(node);
      painted = true;
    }
    return painted;
  }

  function renderMarks() {
    var root = document.querySelector(".main-content");
    if (!root) return;
    var list = marksHere(), missing = 0;
    list.forEach(function (mk) {
      if (root.querySelector('mark[data-id="' + mk.id + '"]')) return;
      var flat = flatten(root);
      var pos = locate(flat, mk);
      if (!pos || !paint(flat, pos.start, pos.end, mk)) missing++;
    });
    if (missing) console.info("[c2a] " + missing + " highlight(s) no longer match this page's text");
  }

  /* ------------------------------ popovers ------------------------------ */

  function closePop() {
    var p = document.getElementById("c2a-pop");
    if (p) p.remove();
  }

  function popAt(x, y) {
    closePop();
    var el = document.createElement("div");
    el.className = "c2a-pop c2a-ui";
    el.id = "c2a-pop";
    document.body.appendChild(el);
    var w = 300;
    el.style.left = Math.max(8, Math.min(x - w / 2, window.innerWidth - w - 12)) + "px";
    el.style.top = (y + window.scrollY + 8) + "px";
    return el;
  }

  function btn(label, cls, fn) {
    var b = document.createElement("button");
    b.textContent = label;
    if (cls) b.className = cls;
    b.addEventListener("click", fn);
    return b;
  }

  /* --------------------------- creating a mark -------------------------- */

  function captureSelection() {
    var sel = window.getSelection();
    if (!sel || sel.isCollapsed || !sel.rangeCount) return null;
    var root = document.querySelector(".main-content");
    var range = sel.getRangeAt(0);
    if (!root || !root.contains(range.commonAncestorContainer)) return null;
    var exact = sel.toString();
    if (!exact.trim() || exact.length < 2) return null;

    var flat = flatten(root);
    var start = offsetOf(flat, range.startContainer, range.startOffset);
    if (start < 0 || flat.text.substr(start, exact.length) !== exact) {
      start = flat.text.indexOf(exact); // fall back for element-boundary selections
    }
    if (start < 0) return null;
    return {
      exact: exact,
      prefix: flat.text.slice(Math.max(0, start - CTX), start),
      suffix: flat.text.slice(start + exact.length, start + exact.length + CTX),
      rect: range.getBoundingClientRect()
    };
  }

  function commitMark(cap, type, comment) {
    var mk = {
      id: uid(), type: type, exact: cap.exact, prefix: cap.prefix, suffix: cap.suffix,
      comment: comment || "", at: new Date().toISOString(), title: pageTitle()
    };
    var list = marksHere();
    list.push(mk);
    setMarksHere(list);
    window.getSelection().removeAllRanges();
    closePop();
    renderMarks();
  }

  function offerSelection() {
    var cap = captureSelection();
    if (!cap) return;
    var el = popAt(cap.rect.left + cap.rect.width / 2, cap.rect.bottom);
    el.appendChild(btn("Highlight", "", function () { commitMark(cap, "note", ""); }));
    el.appendChild(btn("Note…", "", function () { askComment(cap, "note"); }));
    el.appendChild(btn("Flag error…", "danger", function () { askComment(cap, "error"); }));
  }

  function askComment(cap, type) {
    var el = popAt(cap.rect.left + cap.rect.width / 2, cap.rect.bottom);
    var q = document.createElement("div");
    q.className = "quote";
    q.textContent = cap.exact.length > 220 ? cap.exact.slice(0, 220) + "…" : cap.exact;
    el.appendChild(q);
    var ta = document.createElement("textarea");
    ta.placeholder = type === "error"
      ? "What is wrong? (a wrong sign, a bad citation, a broken derivation…)"
      : "Your note";
    el.appendChild(ta);
    el.appendChild(btn("Save", "", function () { commitMark(cap, type, ta.value.trim()); }));
    el.appendChild(btn("Cancel", "", closePop));
    ta.focus();
  }

  /* -------------------------- editing a mark ---------------------------- */

  function issueUrl(mk, path) {
    var page = location.origin + (path || PATH);
    var title = "Erratum: " + (mk.title || "course") + " — " + mk.exact.slice(0, 60).replace(/\s+/g, " ");
    var body = [
      "**Page:** " + page + (mk.id ? "#c2a-" + mk.id : ""),
      "",
      "**Quoted text**",
      "> " + mk.exact.replace(/\n/g, "\n> "),
      "",
      "**What looks wrong**",
      mk.comment || "_(no comment given)_",
      "",
      "---",
      "Reported from the course site."
    ].join("\n");
    return "https://github.com/" + REPO + "/issues/new?labels=erratum"
      + "&title=" + encodeURIComponent(title)
      + "&body=" + encodeURIComponent(body);
  }

  function openMark(id, x, y) {
    var list = marksHere();
    var mk = list.filter(function (m) { return m.id === id; })[0];
    if (!mk) return;
    var el = popAt(x, y);

    if (mk.comment) {
      var c = document.createElement("div");
      c.className = "quote";
      c.textContent = mk.comment;
      el.appendChild(c);
    }
    el.appendChild(btn(mk.type === "error" ? "Unflag" : "Flag as error", "", function () {
      mk.type = mk.type === "error" ? "note" : "error";
      setMarksHere(list);
      document.querySelectorAll('mark[data-id="' + id + '"]').forEach(function (n) {
        n.className = "c2a-mark c2a-" + mk.type;
      });
      closePop();
    }));
    el.appendChild(btn(mk.comment ? "Edit note" : "Add note", "", function () {
      var el2 = popAt(x, y);
      var ta = document.createElement("textarea");
      ta.value = mk.comment || "";
      el2.appendChild(ta);
      el2.appendChild(btn("Save", "", function () {
        mk.comment = ta.value.trim();
        setMarksHere(list);
        closePop();
      }));
      el2.appendChild(btn("Cancel", "", closePop));
      ta.focus();
    }));
    if (mk.type === "error") {
      var a = document.createElement("button");
      a.textContent = "Report on GitHub";
      a.addEventListener("click", function () { window.open(issueUrl(mk), "_blank", "noopener"); });
      el.appendChild(a);
    }
    el.appendChild(btn("Delete", "danger", function () {
      setMarksHere(list.filter(function (m) { return m.id !== id; }));
      document.querySelectorAll('mark[data-id="' + id + '"]').forEach(function (n) {
        var parent = n.parentNode;
        while (n.firstChild) parent.insertBefore(n.firstChild, n);
        parent.removeChild(n);
        parent.normalize();
      });
      closePop();
    }));
  }

  /* ------------------------------ progress ------------------------------ */

  function isUnit() {
    var pages = window.C2A_PAGES || [];
    for (var i = 0; i < pages.length; i++) {
      if (pages[i].url === PATH) return !!pages[i].parent || /orientation/.test(PATH);
    }
    return false;
  }

  function unitList() {
    return (window.C2A_PAGES || []).filter(function (p) {
      return p.parent || /orientation/.test(p.url);
    });
  }

  function mountUnitBar() {
    var root = document.querySelector(".main-content");
    if (!root || !isUnit()) return;
    var h1 = root.querySelector("h1");
    var bar = document.createElement("div");
    bar.id = "c2a-unitbar";
    bar.className = "c2a-ui";
    var b = document.createElement("button");
    var hint = document.createElement("span");
    hint.className = "hint";
    bar.appendChild(b);
    bar.appendChild(hint);

    function sync() {
      var done = !!(getProgress()[PATH] || {}).done;
      bar.classList.toggle("done", done);
      b.textContent = done ? "✓ Completed" : "Mark complete";
      var u = unitList(), n = 0;
      var prog = getProgress();
      u.forEach(function (p) { if ((prog[p.url] || {}).done) n++; });
      hint.textContent = n + " of " + u.length + " units complete"
        + (done && (prog[PATH] || {}).at ? " · finished " + new Date(prog[PATH].at).toLocaleDateString() : "");
      paintProgress();
    }
    b.addEventListener("click", function () {
      var prog = getProgress();
      if ((prog[PATH] || {}).done) delete prog[PATH];
      else prog[PATH] = { done: true, at: new Date().toISOString() };
      save(K_PROGRESS, prog);
      sync();
    });
    if (h1 && h1.nextSibling) h1.parentNode.insertBefore(bar, h1.nextSibling);
    else root.insertBefore(bar, root.firstChild);
    sync();
  }

  function paintProgress() {
    var prog = getProgress();
    document.querySelectorAll(".nav-list-link").forEach(function (a) {
      var href = a.getAttribute("href");
      a.classList.toggle("c2a-done", !!(prog[href] || {}).done);
    });
    var box = document.getElementById("c2a-progress");
    if (!box) return;
    var u = unitList(), n = 0;
    u.forEach(function (p) { if ((prog[p.url] || {}).done) n++; });
    var pct = u.length ? Math.round(100 * n / u.length) : 0;
    box.querySelector("span").textContent = n + " / " + u.length + " units · " + pct + "%";
    box.querySelector("i").style.width = pct + "%";
  }

  function mountProgressMeter() {
    var nav = document.querySelector(".side-bar");
    if (!nav) return;
    var box = document.createElement("div");
    box.id = "c2a-progress";
    box.className = "c2a-ui";
    box.innerHTML = '<span></span><div class="bar"><i></i></div>';
    nav.appendChild(box);
  }

  /* ------------------------------- sidebar ------------------------------ */

  function mountNav() {
    var t = document.createElement("button");
    t.id = "c2a-nav-toggle";
    t.className = "c2a-ui";
    t.title = "Show/hide the navigation panel";
    function sync() {
      var hidden = document.body.classList.contains("c2a-nav-hidden");
      t.textContent = hidden ? "☰ Menu" : "◀ Hide";
    }
    if (load(K_NAVHIDDEN, false)) document.body.classList.add("c2a-nav-hidden");
    t.addEventListener("click", function () {
      var hidden = document.body.classList.toggle("c2a-nav-hidden");
      save(K_NAVHIDDEN, hidden);
      sync();
    });
    document.body.appendChild(t);
    sync();

    // Remember which sections are expanded, across pages.
    var state = load(K_NAV, {});
    document.querySelectorAll(".nav-list-item").forEach(function (li) {
      var exp = li.querySelector(":scope > .nav-list-expander");
      var link = li.querySelector(":scope > .nav-list-link");
      if (!exp || !link) return;
      var key = link.textContent.trim();
      if (Object.prototype.hasOwnProperty.call(state, key)) {
        li.classList.toggle("active", !!state[key]);
      }
      exp.addEventListener("click", function () {
        setTimeout(function () {
          state[key] = li.classList.contains("active");
          save(K_NAV, state);
        }, 0);
      });
    });
  }

  /* ------------------------------- notes -------------------------------- */

  function renderNotes(host) {
    var all = getMarks(), prog = getProgress();
    var paths = Object.keys(all).filter(function (p) { return all[p].length; });
    var total = 0;
    paths.forEach(function (p) { total += all[p].length; });

    var u = unitList(), done = 0;
    u.forEach(function (p) { if ((prog[p.url] || {}).done) done++; });

    host.innerHTML = "";

    var summary = document.createElement("p");
    summary.innerHTML = "<strong>" + done + " of " + u.length + "</strong> units complete · <strong>"
      + total + "</strong> highlight" + (total === 1 ? "" : "s") + " across <strong>"
      + paths.length + "</strong> page" + (paths.length === 1 ? "" : "s") + ".";
    host.appendChild(summary);

    var tb = document.createElement("div");
    tb.className = "c2a-notes-toolbar c2a-ui";
    tb.appendChild(btn("Export JSON", "", function () {
      var blob = new Blob([JSON.stringify({ progress: prog, marks: all }, null, 2)], { type: "application/json" });
      var a = document.createElement("a");
      a.href = URL.createObjectURL(blob);
      a.download = "circuits-to-algorithms-notes.json";
      a.click();
    }));
    tb.appendChild(btn("Import JSON", "", function () {
      var inp = document.createElement("input");
      inp.type = "file"; inp.accept = "application/json";
      inp.addEventListener("change", function () {
        var f = inp.files[0]; if (!f) return;
        var r = new FileReader();
        r.onload = function () {
          try {
            var d = JSON.parse(r.result);
            if (d.marks) save(K_MARKS, d.marks);
            if (d.progress) save(K_PROGRESS, d.progress);
            renderNotes(host);
          } catch (e) { alert("Could not read that file."); }
        };
        r.readAsText(f);
      });
      inp.click();
    }));
    tb.appendChild(btn("Clear everything", "danger", function () {
      if (confirm("Delete all highlights and progress in this browser? This cannot be undone.")) {
        localStorage.removeItem(K_MARKS);
        localStorage.removeItem(K_PROGRESS);
        renderNotes(host);
      }
    }));
    host.appendChild(tb);

    if (!total) {
      var p = document.createElement("p");
      p.className = "c2a-empty";
      p.textContent = "No highlights yet. Select any text in a unit and choose Highlight, Note, or Flag error.";
      host.appendChild(p);
      return;
    }

    paths.sort();
    paths.forEach(function (path) {
      var g = document.createElement("div");
      g.className = "c2a-group";
      var h = document.createElement("h3");
      var a = document.createElement("a");
      a.href = path;
      a.textContent = all[path][0].title || path;
      h.appendChild(a);
      g.appendChild(h);

      all[path].forEach(function (mk) {
        var e = document.createElement("div");
        e.className = "c2a-entry " + mk.type;
        var q = document.createElement("blockquote");
        q.textContent = mk.exact;
        e.appendChild(q);
        if (mk.comment) {
          var c = document.createElement("p");
          c.className = "comment";
          c.textContent = mk.comment;
          e.appendChild(c);
        }
        var meta = document.createElement("div");
        meta.className = "meta";
        var badge = document.createElement("span");
        badge.className = "c2a-badge " + mk.type;
        badge.textContent = mk.type === "error" ? "error" : "note";
        meta.appendChild(badge);
        var when = document.createElement("span");
        when.textContent = new Date(mk.at).toLocaleDateString();
        meta.appendChild(when);
        var go = document.createElement("a");
        go.href = path + "#c2a-" + mk.id;
        go.textContent = "Go to passage";
        meta.appendChild(go);
        if (mk.type === "error") {
          var rep = document.createElement("a");
          rep.href = issueUrl(mk, path);
          rep.target = "_blank";
          rep.rel = "noopener";
          rep.textContent = "Report on GitHub →";
          meta.appendChild(rep);
        }
        meta.appendChild(btn("Delete", "", function () {
          var d = getMarks();
          d[path] = (d[path] || []).filter(function (m) { return m.id !== mk.id; });
          if (!d[path].length) delete d[path];
          save(K_MARKS, d);
          renderNotes(host);
        }));
        e.appendChild(meta);
        g.appendChild(e);
      });
      host.appendChild(g);
    });
  }

  /* -------------------------------- boot -------------------------------- */

  function boot() {
    mountNav();
    mountProgressMeter();
    mountUnitBar();
    paintProgress();
    renderMarks();

    var host = document.getElementById("c2a-notes");
    if (host) renderNotes(host);

    document.addEventListener("mouseup", function (ev) {
      if (ev.target.closest && ev.target.closest(".c2a-pop")) return;
      setTimeout(function () {
        var m = ev.target.closest && ev.target.closest("mark.c2a-mark");
        if (m && window.getSelection().isCollapsed) {
          openMark(m.dataset.id, ev.clientX, ev.clientY);
        } else if (!window.getSelection().isCollapsed) {
          offerSelection();
        } else {
          closePop();
        }
      }, 10);
    });
    document.addEventListener("keydown", function (e) { if (e.key === "Escape") closePop(); });

    if (location.hash.indexOf("#c2a-") === 0) {
      var el = document.getElementById(location.hash.slice(1));
      if (el) { el.scrollIntoView({ block: "center" }); el.classList.add("c2a-flash"); }
    }
  }

  // Wait for MathJax: it rewrites the DOM, and anchors must be computed
  // against the post-typeset text or offsets will not match on reload.
  function ready() {
    if (window.MathJax && window.MathJax.startup && window.MathJax.startup.promise) {
      window.MathJax.startup.promise.then(boot).catch(boot);
    } else {
      boot();
    }
  }
  if (document.readyState === "complete") ready();
  else window.addEventListener("load", ready);

  window.C2A = { renderNotes: renderNotes, issueUrl: issueUrl };
})();
