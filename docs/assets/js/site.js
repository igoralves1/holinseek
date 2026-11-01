/* HolinSeek — comportamento compartilhado do site.
 *
 *  1. tema claro/escuro, lembrado no navegador;
 *  2. marcação da página atual na navegação;
 *  3. desenho das estruturas moleculares (SmilesDrawer).
 *
 * MathJax é configurado no próprio HTML, antes de a biblioteca carregar.
 */

(function () {
  "use strict";

  /* ---------------- tema ---------------- */

  var KEY = "holinseek-theme";

  function applyTheme(theme) {
    if (theme === "dark" || theme === "light") {
      document.documentElement.setAttribute("data-theme", theme);
    } else {
      document.documentElement.removeAttribute("data-theme");
    }
  }

  function storedTheme() {
    try {
      return localStorage.getItem(KEY);
    } catch (e) {
      return null;
    }
  }

  function storeTheme(theme) {
    try {
      if (theme) localStorage.setItem(KEY, theme);
      else localStorage.removeItem(KEY);
    } catch (e) {
      /* navegação privada ou armazenamento bloqueado: segue sem lembrar */
    }
  }

  applyTheme(storedTheme());

  function currentTheme() {
    var explicit = document.documentElement.getAttribute("data-theme");
    if (explicit) return explicit;
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches
      ? "dark"
      : "light";
  }

  function setupThemeToggle() {
    var btn = document.querySelector(".theme-toggle");
    if (!btn) return;

    function label() {
      btn.textContent = currentTheme() === "dark" ? "claro" : "escuro";
      btn.setAttribute(
        "aria-label",
        currentTheme() === "dark" ? "Mudar para tema claro" : "Mudar para tema escuro"
      );
    }

    label();
    btn.addEventListener("click", function () {
      var next = currentTheme() === "dark" ? "light" : "dark";
      applyTheme(next);
      storeTheme(next);
      label();
      drawMolecules();
    });
  }

  /* ---------------- navegação ---------------- */

  function markCurrentPage() {
    var here = location.pathname.split("/").pop() || "index.html";
    document.querySelectorAll("nav.main a").forEach(function (a) {
      var target = a.getAttribute("href");
      if (target === here || (here === "" && target === "index.html")) {
        a.setAttribute("aria-current", "page");
      }
    });
  }

  /* ---------------- estruturas moleculares ---------------- */

  /* Cada <canvas data-smiles="..."> vira um desenho 2D da molécula.
   * Se a biblioteca não carregar (sem rede, CDN bloqueado) ou o SMILES não
   * for aceito, o canvas some e o texto de reserva aparece — a página
   * continua legível, porque a fórmula molecular está sempre escrita embaixo.
   */
  var molSeq = 0;

  function drawMolecules() {
    var canvases = document.querySelectorAll("canvas[data-smiles]");
    if (!canvases.length) return;

    if (typeof window.SmilesDrawer === "undefined" || !window.SmilesDrawer.SmiDrawer) {
      canvases.forEach(showFallback);
      return;
    }

    /* As moléculas são sempre desenhadas sobre fundo claro: o cartão que as
     * contém tem fundo branco nos dois temas, porque estrutura química com
     * traço claro sobre escuro fica difícil de ler. */
    canvases.forEach(function (canvas) {
      if (!canvas.id) canvas.id = "mol-" + ++molSeq;
      canvas.style.display = "block";

      var drawer = new window.SmilesDrawer.SmiDrawer({
        width: canvas.width,
        height: canvas.height,
        bondThickness: 1.1,
        padding: 14,
        terminalCarbons: false
      });

      try {
        drawer.draw(canvas.getAttribute("data-smiles"), "#" + canvas.id, "light");
      } catch (err) {
        showFallback(canvas);
        return;
      }

      /* O desenho passa por uma conversão assíncrona de SVG para canvas, então
       * só dá para saber se deu certo olhando os pixels um instante depois. */
      window.setTimeout(function () {
        if (isBlank(canvas)) showFallback(canvas);
      }, 1200);
    });
  }

  function isBlank(canvas) {
    try {
      var data = canvas.getContext("2d").getImageData(0, 0, canvas.width, canvas.height).data;
      var painted = 0;
      for (var i = 3; i < data.length; i += 4) {
        if (data[i] > 0 && ++painted > 150) return false;
      }
      return true;
    } catch (e) {
      return false; /* na dúvida, mantém o desenho */
    }
  }

  function showFallback(canvas) {
    canvas.style.display = "none";
    var fallback = canvas.parentElement
      ? canvas.parentElement.querySelector(".mol-fallback")
      : null;
    if (fallback) fallback.style.display = "block";
  }

  /* ---------------- início ---------------- */

  function init() {
    setupThemeToggle();
    markCurrentPage();
    drawMolecules();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }

  window.addEventListener("load", drawMolecules);
})();
