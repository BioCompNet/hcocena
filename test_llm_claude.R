# =====================================================================
# Test: hc_module_function_llm() mit eigener Gen-Liste über Claude
# ---------------------------------------------------------------------
# Aufruf:  Rscript test_llm_claude.R
# (oder Zeile fuer Zeile in RStudio / R-Konsole ausfuehren)
# =====================================================================

# --- Paket laden (Entwicklungsmodus, aus dem Repo-Root) --------------
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(hcocena)
}

# --- 1) DEIN API-Key -------------------------------------------------
# Variante A: vorher in der Shell setzen, z. B.
#   PowerShell:  $env:ANTHROPIC_API_KEY = "sk-ant-..."
# Variante B: hier direkt eintragen (NICHT eingecheckt lassen!):
# Sys.setenv(ANTHROPIC_API_KEY = "sk-ant-...")

api_key <- Sys.getenv("ANTHROPIC_API_KEY")
if (!nzchar(api_key)) {
  stop("Kein ANTHROPIC_API_KEY gesetzt. Bitte oben Variante A oder B nutzen.")
}

# --- 2) DEINE Gen-Liste (hier ersetzen) ------------------------------
my_genes <- c(
  "STAT1", "IRF7", "CXCL10", "GBP1", "IFI44L", "ISG15", "MX1", "OAS1"
)

# Optionaler biologischer Kontext (steuert die Interpretation)
my_context <- "Interferon-getriebenes Blutmodul bei akuter Virusinfektion"

# --- 3) Lauf gegen Claude --------------------------------------------
res <- hc_module_function_llm(
  genes       = my_genes,
  context     = my_context,        # optional, kann auch weggelassen werden
  label       = "mein_test",       # frei waehlbar
  llm         = "claude",
  api_key     = api_key,
  # claude_model = "claude-opus-4-8",  # optional; Default ist claude-sonnet-4-6
  save_to_hc  = FALSE,             # wichtig: kein hc-Objekt noetig
  verbose     = TRUE
)

# --- 4) Ergebnis ausgeben --------------------------------------------
cat("\n================ ERGEBNIS ================\n")
cat("Provider/Model :", res$llm, "/", res$model, "\n")
cat("Gene gesendet  :", res$gene_count_sent, "von", res$gene_count_input, "\n")
cat("Status         :", res$status, "\n\n")

resp <- res$response
cat("general_processes:\n  ", resp$general_processes, "\n\n")
cat("contextual_state :\n  ", resp$contextual_state, "\n\n")
cat("key_regulators   :\n  ", resp$key_regulators, "\n")
cat("==========================================\n")

# Roh-Antwort fuers Debugging:
# cat(res$raw_response_text, "\n")
