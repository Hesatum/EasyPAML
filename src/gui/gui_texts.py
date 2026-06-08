"""
gui_texts.py — Dicionário Central de Textos da Interface EasyPAML
=================================================================

Internacionalização (i18n):
  - TEXTS_PT  : textos em Português do Brasil (padrão)
  - TEXTS_EN  : textos em Inglês
  - TEXTS     : proxy transparente — roteia para o idioma ativo.
                Todo código existente que usa TEXTS["key"] continua
                funcionando sem nenhuma alteração.

Para mudar o idioma em tempo de execução:
    from src.gui.gui_texts import set_language
    set_language('en')   # ou 'pt'

Para ler o idioma atual:
    from src.gui.gui_texts import get_language
    get_language()  # → 'pt' ou 'en'

ORGANIZAÇÃO DOS DICIONÁRIOS
----------------------------
  Seção 1 — main_gui.py › ModelConfigWindow
  Seção 2 — main_gui.py › TreeLabelWindow
  Seção 3 — main_gui.py › App  (janela principal, sidebar, log)
  Seção 4 — results_viewer.py › ResultsViewerWindow

TEMPLATES (strings com {placeholders})
---------------------------------------
  Use .format() para substituir variáveis dinâmicas. Exemplos:
    TEXTS["sites_file_not_found"].format(filename="gene_M8.txt")
    TEXTS["lrt_footer_template"].format(total=30, sig=5, df=2)
    TEXTS["status_stops_template"].format(n=3)
"""

# ═══════════════════════════════════════════════════════════════════════════
# PORTUGUÊS DO BRASIL
# ═══════════════════════════════════════════════════════════════════════════

TEXTS_PT: dict[str, object] = {

    # ── Seção 1 — ModelConfigWindow ─────────────────────────────────────────
    "model_config_header":      "Modelo: {model_code}",
    "model_config_btn_cancel":  "Cancelar",
    "model_config_btn_save":    "Salvar",
    "model_status_default":     "Default",
    "model_status_configured":  "Configurado",

    # ── Seção 2 — TreeLabelWindow ────────────────────────────────────────────
    "tree_labeler_sidebar_title": "Instrucoes",

    "tree_labeler_instructions_branchsite": (
        "• Clique nos círculos\n  para marcar/desmarcar\n  foreground (#1).\n\n"
        "• O modo branch-site só permite uma tag #1\n.\n\n"
        "• Cor: Vermelho"
    ),
    "tree_labeler_instructions_branch": (
        "• Clique nos círculos\n  para atribuir tags e testar o ômega dos ramos.\n\n"
        "• Digite número\n  da tag (1, 2, 3...).\n\n"
        "• Cada número\n  recebe cor única."
    ),

    "tree_labeler_legend_title": "Tags Ativas",
    "tree_labeler_no_tags":      "Nenhuma tag ativa",
    "tree_labeler_btn_save":     "Salvar",
    "tree_labeler_btn_cancel":   "Cancelar",

    # ── Seção 3 — App (janela principal) ────────────────────────────────────
    "app_sidebar_title":    "EasyPAML",
    "app_sidebar_subtitle": "Seleção Positiva",

    "section_files":   "ARQUIVOS",
    "section_results": "RESULTADOS",
    "section_config":  "CONFIGURACOES",

    "btn_input_folder":  "Pasta de alinhamentos",
    "btn_tree_file":     "Tree file",
    "btn_output_folder": "Pasta de Resultados",
    "label_not_selected": "Nao selecionado",

    "btn_view_results":          "Ver Resultados",
    "btn_update_results":        "Atualizar Resultados",
    "label_update_results_hint": "Atualizar arquivos de analise",

    "label_omega_initial":    "dN/dS Inicial (w):",
    "label_remove_gaps":      "Limpar alinhamento antes de analisar",
    "label_remove_gaps_hint": (
        "Recomendado para a maioria das análises.\n"
        "Remove automaticamente colunas do alinhamento que\n"
        "contenham gaps (—), bases incertas (N, ?) ou códons\n"
        "incompletos — dados assim podem distorcer as\n"
        "estimativas de seleção positiva.\n"
        "Desative apenas se o alinhamento já foi curado\n"
        "manualmente e você não quer perder nenhum sítio."
    ),
    "label_cpus": "CPUs (paralelismo):",

    "label_ignore_stops":      "Ignorar Stop Codons",
    "label_ignore_stops_hint": (
        "Envia Enter automaticamente ao encontrar stop codons.\n"
        "No modo paralelo (CPUs > 1) é sempre ativado."
    ),

    "label_auto_prune":      "Poda Automática de Árvore",
    "label_auto_prune_hint": (
        "Poda a árvore para incluir apenas os\n"
        "taxons presentes em cada gene automaticamente."
    ),

    "tab_site_models":  "Site Models",
    "tab_branch_model": "Branch Model",
    "tab_branchsite":   "Branch-Site",

    "status_ready":   "● Pronto",
    "status_running": "● Executando",
    "status_paused":  "● Pausada",
    "status_stopped": "● Parada",
    "status_stops_template": "■ Stops: {n}",
    "label_neutral_models":  "Modelos nulos automáticos",

    "btn_run":           "Iniciar",
    "btn_pause":         "Pausar",
    "btn_resume":        "Retomar",
    "btn_stop":          "Parar",
    "btn_resolve_stops": "Ignorar Stop codons",

    "log_header_title":    "LOG DE EXECUCAO",
    "log_header_subtitle": "output em tempo real",

    "log_welcome": (
        "  ╔══════════════════════════════════════════════╗\n"
        "  ║   EasyPAML  ·  PAML / CODEML interface      ║\n"
        "  ╚══════════════════════════════════════════════╝\n\n"
        "  1. Selecione a pasta de alinhamentos (.fas / .fasta / .phy / .phylip)\n"
        "  2. Escolha o Tree File (.tree / .tre / .nwk / .txt)\n"
        "  3. Defina a pasta de saída de resultados\n"
        "  4. Marque os modelos e clique  ▶ Iniciar\n\n"
        "  ─────────────────────────────────────────────────\n\n"
    ),

    "btn_label_branch":     "Marcar Branch (Multiplos Ramos)",
    "btn_label_branchsite": "Marcar Branch-site",

    "model_info_test_type":      "Tipo de Teste:",
    "model_info_params":         "Parametros:",
    "model_info_purpose":        "Proposito:",
    "model_info_interpretation": "Interpretacao:",
    "model_info_use_case":       "Quando usar:",
    "model_info_references":     "Referencias:",

    # ── Seção 4 — ResultsViewerWindow ───────────────────────────────────────
    "viewer_window_title":       "EasyPAML - Painel de Analise",
    "viewer_error_no_tsv":       "Arquivo analysis_summary.tsv não encontrado!",
    "viewer_error_run_analysis": "Execute uma análise para gerar resultados.",

    "viewer_header_title":    "EasyPAML  —  Resultados",
    "viewer_header_subtitle": "Análise de Seleção Positiva  ·  CODEML / PAML",

    "viewer_tab_lrt":                "LRT & P-values",
    "viewer_tab_omega":              "ω > 1 Global",
    "viewer_tab_sites":              "Sítios Positivos",
    "viewer_tab_branchsite_classes": "Branch-site Classes",
    "viewer_tab_branch":             "Branch Analysis",
    "viewer_tab_export":             "Exportar",

    "stats_total_genes":        "Total de Genes",
    "stats_models_run":         "Modelos Rodados",
    "stats_positive_selection": "Seleção Positiva Global",
    "stats_avg_omega":          "ω Médio",

    "lrt_no_comparisons": "[!] Sem comparacoes LRT disponiveis",
    "lrt_label_model":    "Modelo testado:",

    "lrt_branch_warning": (
        "Branch e Branch-site requerem sequencias com mais de "
        "200 pb para estimativas confiaveis de ω.  Genes mais curtos podem "
        "produzir estimativas instaveis ou nao convergir."
    ),

    "lrt_no_data_for_comparison": "[!] Nenhum dado LRT disponivel para esta comparacao",

    "lrt_footer_template": (
        "Total: {total} genes  ·  "
        "Significantes p < 0.05: {sig}  ·  "
        "df = {df}  ·  "
        "* = significante  ·  * ω>1 = selecao positiva"
    ),

    "pos_sel_tab_title": "Selecao Global — ω > 1 no gene inteiro",
    "pos_sel_tab_criterion": (
        "Critério: ω médio do modelo M2a ou M8 > 1.0  AND  LRT p < 0.05  ·  "
        "Diferente de seleção em sítios específicos (aba Sítios Positivos)"
    ),
    "pos_sel_none_found":      "Nenhum sinal de seleção positiva detectado",
    "pos_sel_criterion_short": "Critério: ω > 1.0  AND  p-valor < 0.05",
    "pos_sel_badge":           "* Positivo",

    "sites_label_model":    "Modelo:",
    "sites_label_gene":     "Gene:",
    "sites_label_analysis": "Análise:",
    "sites_label_filter":   "Filtrar Pr(w>1) ≥",
    "sites_btn_update":     "Atualizar",

    "sites_file_not_found": "[X] Arquivo nao encontrado: {filename}",
    "sites_parse_error":    "[Erro] Erro ao parsear arquivo:\n{error}",
    "sites_no_sites":       "Nenhum sitio acima do limiar Pr(w>1) >= {threshold}",

    "sites_table_headers": [
        "Posição",    # 70 px
        "AA",         # 55 px
        "Pr(ω>1)",    # 100 px
        "Sig",        # 55 px
        "ω (média)",  # 105 px
        "ω − SE",     # 105 px
        "ω + SE",     # 105 px
    ],

    "branchsite_classes_gene_label": "Gene:",
    "branchsite_classes_not_found":  "Gene não encontrado",
    "branchsite_classes_header":     "Classes do Modelo Branch-site — {gene}",

    "branchsite_classes_table_headers": [
        "Classe",        # 120 px
        "Proporção",     # 150 px
        "Background ω",  # 150 px
        "Foreground ω",  # 150 px
    ],

    "branchsite_classes_bg_w_placeholder": "Ver arquivo",
    "branchsite_classes_footer": (
        "Valores de Foreground ω: valores elevados (> 1.0) indicam seleção positiva "
        "no ramo foreground para aquela classe de sítios"
    ),

    "branch_tab_title":  "Branch Analysis — Cladograma por dN/dS",
    "branch_tab_legend": (
        "Vermelho (w<1) · Amarelo (w=1) · Azul (w>1)"
        "   |   Clique num nó interno para girar"
    ),
    "branch_no_data_title": "Sem dados do Branch Model",
    "branch_no_data_hint":  "Execute o modelo Branch com uma arvore etiquetada.",
    "branch_label_gene":     "Gene:",
    "branch_label_outgroup": "Outgroup:",
    "branch_outgroup_none":  "(nenhum)",
    "branch_btn_export_png": "Exportar PNG",

    "export_tab_title": "Exportar Resultados",
    "export_btn":       "Exportar →",

    "export_options": [
        ("Excel (.xlsx)",    "Tabela com formatação profissional"),
        ("CSV",              "Formato universal compatível"),
        ("Gráficos (PNG)",   "Exportar gráficos em alta resolução"),
        ("Relatório (HTML)", "Relatório completo interativo"),
    ],

    # ── TreeLabelWindow — mensagens de erro inline ───────────────────────
    "tree_err_no_biopython":  "[!] Biopython não instalado. Execute: pip install biopython",
    "tree_err_no_tree":       "[!] Nenhuma arvore selecionada.",
    "tree_err_load":          "[Erro] Erro ao carregar arvore:\n{error}",

    # ── TreeLabelWindow — diálogos de tag ────────────────────────────────
    "tag_dialog_edit_title":  "Editar Tag",
    "tag_dialog_edit_prompt": "Ramo atual: {tag}\n\nDigite novo número ou 'remover':",
    "tag_dialog_new_title":   "Número da Tag",
    "tag_dialog_new_prompt":  "Digite o número da tag:\n(ex: 1 para #1, 2 para #2)",

    # ── App — CPU label ──────────────────────────────────────────────────
    "label_cpus_detected":    "(detectado: {n} núcleos)",

    # ── App — mensagens de log ───────────────────────────────────────────
    "log_no_tree_selected":   "[!] Selecione uma arvore (.nwk) primeiro.\n",
    "log_no_output_folder":   "[!] Selecione uma pasta de saida primeiro.\n",
    "log_updating_results":   ">> ATUALIZANDO RESULTADOS\n",
    "log_detecting_models":   "Detectando modelos presentes...\n",
    "log_lrt_comparisons":    "Determinando comparações para LRT:\n",
    "log_lrt_M0_M1a":         "  • M0 (null) vs M1a (alt) - Variação de ω entre sítios\n",
    "log_lrt_M1a_M2a":        "  • M1a (null) vs M2a (alt) - Seleção positiva\n",
    "log_lrt_M7_M8":          "  • M7 (null) vs M8 (alt) - Seleção positiva (Beta)\n",
    "log_lrt_M0_Branch":      "  • M0 (null) vs Branch (alt) - Seleção por ramo\n",
    "log_lrt_BranchSite":     "  • Branch-site_null (null) vs Branch-site (alt) - Seleção branch-site\n",
    "log_lrt_total":          "Total de {n} comparações encontradas.\n\n",
    "log_regenerating":       "Regenerando arquivos de síntese...\n",
    "log_update_done":        "\n[OK] ATUALIZACAO CONCLUIDA COM SUCESSO!\n",
    "log_resolving_stops":    "▶ Resolvendo todos os stops deste Exon automaticamente...\n",
    "log_analysis_start":     "▶ INICIANDO ANALISE BATCH\n",
    "log_analysis_stopped":   "■ ANALISE INTERROMPIDA\n",
    "log_analysis_done":      "[OK] ANALISE CONCLUIDA\n",

    # ── App — diálogos de arquivo ────────────────────────────────────────
    "dialog_select_results_folder": "Selecione a pasta com resultados para atualizar síntese",

    # ── App — troca de idioma ────────────────────────────────────────────
    "lang_switch_message":  "Reinicie o EasyPAML para aplicar o novo idioma.",
    "lang_switch_confirm":  "Reiniciar agora?",
    "lang_switch_err":      "Não foi possível reiniciar automaticamente:\n{error}\n\nReabra manualmente.",

    # ── ResultsViewerWindow — mensagens inline ───────────────────────────
    "viewer_genes_loaded":      "{n} genes carregados",
    "viewer_gene_not_found":    "Gene nao encontrado.",
    "viewer_branch_no_file":    "Arquivo de resultados Branch nao encontrado.",
    "viewer_branch_read_err":   "Erro ao ler tabela: {error}",
    "viewer_branch_no_table":   "Tabela dN & dS nao encontrada no arquivo.",
    "viewer_branch_invalid":    "Dados invalidos.",
    "viewer_lrt_parse_err":     "[Erro] Erro ao parsear comparacao",
    "viewer_sites_subtitle":    "Modelo: {model}  ·  Análise: {method}  ·  {omega}",
    "viewer_sites_count":       "{n} sítios",

    # ── ResultsViewerWindow — caixas de mensagem ─────────────────────────
    "msg_warning":              "Aviso",
    "msg_success":              "Sucesso",
    "msg_error":                "Erro",
    "msg_no_figure":            "Nenhuma figura para exportar.",
    "msg_exported_to":          "Exportado para:\n{path}",
    "msg_export_err":           "Erro ao exportar:\n{error}",
    "msg_no_lrt":               "Nenhum resultado LRT encontrado.",
    "msg_excel_exported":       "Excel exportado com {n} aba(s):\n{path}",
    "msg_excel_err":            "Erro ao exportar Excel:\n{error}",
    "msg_csv_exported":         "{n} arquivo(s) exportado(s):\n{files}",
    "msg_no_csv_data":          "Nenhum dado disponível para exportar.",
    "msg_csv_err":              "Erro ao exportar CSV:\n{error}",
    "msg_html_exported":        "Relatório HTML exportado:\n{path}",
    "msg_html_err":             "Erro ao exportar HTML:\n{error}",

    # ── ResultsViewerWindow — diálogos de arquivo ────────────────────────
    "dialog_export_cladogram":  "Exportar cladograma",
    "dialog_save_csv":          "Salvar CSVs — escolha o nome base (sem extensão)",

    # ── Gráficos (matplotlib) ─────────────────────────────────────────────
    "chart_omega_dist":         "Distribuição de ω",
    "chart_freq":               "Frequência",
}


# ═══════════════════════════════════════════════════════════════════════════
# ENGLISH
# ═══════════════════════════════════════════════════════════════════════════

TEXTS_EN: dict[str, object] = {

    # ── Section 1 — ModelConfigWindow ───────────────────────────────────────
    "model_config_header":      "Model: {model_code}",
    "model_config_btn_cancel":  "Cancel",
    "model_config_btn_save":    "Save",
    "model_status_default":     "Default",
    "model_status_configured":  "Configured",

    # ── Section 2 — TreeLabelWindow ─────────────────────────────────────────
    "tree_labeler_sidebar_title": "Instructions",

    "tree_labeler_instructions_branchsite": (
        "• Click on circles\n  to mark/unmark\n  foreground (#1).\n\n"
        "• Branch-site mode allows only one #1 tag.\n\n"
        "• Color: Red"
    ),
    "tree_labeler_instructions_branch": (
        "• Click on circles\n  to assign tags and test branch ω.\n\n"
        "• Type the tag number\n  (1, 2, 3...).\n\n"
        "• Each number\n  receives a unique color."
    ),

    "tree_labeler_legend_title": "Active Tags",
    "tree_labeler_no_tags":      "No active tags",
    "tree_labeler_btn_save":     "Save",
    "tree_labeler_btn_cancel":   "Cancel",

    # ── Section 3 — App (main window) ───────────────────────────────────────
    "app_sidebar_title":    "EasyPAML",
    "app_sidebar_subtitle": "Positive Selection",

    "section_files":   "FILES",
    "section_results": "RESULTS",
    "section_config":  "SETTINGS",

    "btn_input_folder":   "Alignments folder",
    "btn_tree_file":      "Tree file",
    "btn_output_folder":  "Output folder",
    "label_not_selected": "Not selected",

    "btn_view_results":          "View Results",
    "btn_update_results":        "Update Results",
    "label_update_results_hint": "Regenerate analysis files",

    "label_omega_initial":    "Initial dN/dS (w):",
    "label_remove_gaps":      "Clean alignment before analysis",
    "label_remove_gaps_hint": (
        "Recommended for most analyses.\n"
        "Automatically removes alignment columns that\n"
        "contain gaps (—), ambiguous bases (N, ?) or\n"
        "incomplete codons — such data can distort\n"
        "positive selection estimates.\n"
        "Disable only if the alignment has already been\n"
        "manually curated and no sites should be removed."
    ),
    "label_cpus": "CPUs (parallelism):",

    "label_ignore_stops":      "Ignore Stop Codons",
    "label_ignore_stops_hint": (
        "Automatically sends Enter when stop codons are found.\n"
        "In parallel mode (CPUs > 1) always enabled."
    ),

    "label_auto_prune":      "Auto-prune Tree",
    "label_auto_prune_hint": (
        "Prunes the tree to include only the\n"
        "taxa present in each gene automatically."
    ),

    "tab_site_models":  "Site Models",
    "tab_branch_model": "Branch Model",
    "tab_branchsite":   "Branch-Site",

    "status_ready":   "● Ready",
    "status_running": "● Running",
    "status_paused":  "● Paused",
    "status_stopped": "● Stopped",
    "status_stops_template": "■ Stops: {n}",
    "label_neutral_models":  "Automatic null models",

    "btn_run":           "Run",
    "btn_pause":         "Pause",
    "btn_resume":        "Resume",
    "btn_stop":          "Stop",
    "btn_resolve_stops": "Ignore Stop Codons",

    "log_header_title":    "EXECUTION LOG",
    "log_header_subtitle": "real-time output",

    "log_welcome": (
        "  ╔══════════════════════════════════════════════╗\n"
        "  ║   EasyPAML  ·  PAML / CODEML interface      ║\n"
        "  ╚══════════════════════════════════════════════╝\n\n"
        "  1. Select the alignments folder (.fas / .fasta / .phy / .phylip)\n"
        "  2. Choose the Tree File (.tree / .tre / .nwk / .txt)\n"
        "  3. Set the output folder\n"
        "  4. Select models and click  ▶ Run\n\n"
        "  ─────────────────────────────────────────────────\n\n"
    ),

    "btn_label_branch":     "Label Branch (Multiple Branches)",
    "btn_label_branchsite": "Label Branch-site",

    "model_info_test_type":      "Test Type:",
    "model_info_params":         "Parameters:",
    "model_info_purpose":        "Purpose:",
    "model_info_interpretation": "Interpretation:",
    "model_info_use_case":       "When to use:",
    "model_info_references":     "References:",

    # ── Section 4 — ResultsViewerWindow ─────────────────────────────────────
    "viewer_window_title":       "EasyPAML - Analysis Panel",
    "viewer_error_no_tsv":       "File analysis_summary.tsv not found!",
    "viewer_error_run_analysis": "Run an analysis to generate results.",

    "viewer_header_title":    "EasyPAML  —  Results",
    "viewer_header_subtitle": "Positive Selection Analysis  ·  CODEML / PAML",

    "viewer_tab_lrt":                "LRT & P-values",
    "viewer_tab_omega":              "ω > 1 Global",
    "viewer_tab_sites":              "Positive Sites",
    "viewer_tab_branchsite_classes": "Branch-site Classes",
    "viewer_tab_branch":             "Branch Analysis",
    "viewer_tab_export":             "Export",

    "stats_total_genes":        "Total Genes",
    "stats_models_run":         "Models Run",
    "stats_positive_selection": "Global Positive Selection",
    "stats_avg_omega":          "Average ω",

    "lrt_no_comparisons": "[!] No LRT comparisons available",
    "lrt_label_model":    "Tested model:",

    "lrt_branch_warning": (
        "Branch and Branch-site require sequences longer than "
        "200 bp for reliable ω estimates. Shorter genes may "
        "produce unstable estimates or fail to converge."
    ),

    "lrt_no_data_for_comparison": "[!] No LRT data available for this comparison",

    "lrt_footer_template": (
        "Total: {total} genes  ·  "
        "Significant p < 0.05: {sig}  ·  "
        "df = {df}  ·  "
        "* = significant  ·  * ω>1 = positive selection"
    ),

    "pos_sel_tab_title": "Global Selection — ω > 1 across the whole gene",
    "pos_sel_tab_criterion": (
        "Criterion: mean ω from M2a or M8 > 1.0  AND  LRT p < 0.05  ·  "
        "Different from site-specific selection (Positive Sites tab)"
    ),
    "pos_sel_none_found":      "No positive selection signal detected",
    "pos_sel_criterion_short": "Criterion: ω > 1.0  AND  p-value < 0.05",
    "pos_sel_badge":           "* Positive",

    "sites_label_model":    "Model:",
    "sites_label_gene":     "Gene:",
    "sites_label_analysis": "Analysis:",
    "sites_label_filter":   "Filter Pr(w>1) ≥",
    "sites_btn_update":     "Update",

    "sites_file_not_found": "[X] File not found: {filename}",
    "sites_parse_error":    "[Error] Error parsing file:\n{error}",
    "sites_no_sites":       "No sites above threshold Pr(w>1) >= {threshold}",

    "sites_table_headers": [
        "Position",  # 70 px
        "AA",        # 55 px
        "Pr(ω>1)",   # 100 px
        "Sig",       # 55 px
        "ω (mean)",  # 105 px
        "ω − SE",    # 105 px
        "ω + SE",    # 105 px
    ],

    "branchsite_classes_gene_label": "Gene:",
    "branchsite_classes_not_found":  "Gene not found",
    "branchsite_classes_header":     "Branch-site Model Classes — {gene}",

    "branchsite_classes_table_headers": [
        "Class",        # 120 px
        "Proportion",   # 150 px
        "Background ω", # 150 px
        "Foreground ω", # 150 px
    ],

    "branchsite_classes_bg_w_placeholder": "See file",
    "branchsite_classes_footer": (
        "Foreground ω values: high values (> 1.0) indicate positive selection "
        "on the foreground branch for that site class"
    ),

    "branch_tab_title":  "Branch Analysis — Cladogram by dN/dS",
    "branch_tab_legend": (
        "Red (w<1) · Yellow (w=1) · Blue (w>1)"
        "   |   Click an internal node to rotate"
    ),
    "branch_no_data_title": "No Branch Model data",
    "branch_no_data_hint":  "Run the Branch model with a labeled tree.",
    "branch_label_gene":     "Gene:",
    "branch_label_outgroup": "Outgroup:",
    "branch_outgroup_none":  "(none)",
    "branch_btn_export_png": "Export PNG",

    "export_tab_title": "Export Results",
    "export_btn":       "Export →",

    "export_options": [
        ("Excel (.xlsx)",   "Table with professional formatting"),
        ("CSV",             "Universal compatible format"),
        ("Charts (PNG)",    "Export charts in high resolution"),
        ("Report (HTML)",   "Complete interactive report"),
    ],

    # ── TreeLabelWindow — inline error messages ──────────────────────────
    "tree_err_no_biopython":  "[!] Biopython not installed. Run: pip install biopython",
    "tree_err_no_tree":       "[!] No tree selected.",
    "tree_err_load":          "[Error] Failed to load tree:\n{error}",

    # ── TreeLabelWindow — tag dialogs ────────────────────────────────────
    "tag_dialog_edit_title":  "Edit Tag",
    "tag_dialog_edit_prompt": "Current branch: {tag}\n\nEnter new number or 'remove':",
    "tag_dialog_new_title":   "Tag Number",
    "tag_dialog_new_prompt":  "Enter tag number:\n(e.g. 1 for #1, 2 for #2)",

    # ── App — CPU label ──────────────────────────────────────────────────
    "label_cpus_detected":    "(detected: {n} cores)",

    # ── App — log messages ───────────────────────────────────────────────
    "log_no_tree_selected":   "[!] Please select a tree file (.nwk) first.\n",
    "log_no_output_folder":   "[!] Please select an output folder first.\n",
    "log_updating_results":   ">> UPDATING RESULTS\n",
    "log_detecting_models":   "Detecting present models...\n",
    "log_lrt_comparisons":    "Determining LRT comparisons:\n",
    "log_lrt_M0_M1a":         "  • M0 (null) vs M1a (alt) - ω variation across sites\n",
    "log_lrt_M1a_M2a":        "  • M1a (null) vs M2a (alt) - Positive selection\n",
    "log_lrt_M7_M8":          "  • M7 (null) vs M8 (alt) - Positive selection (Beta)\n",
    "log_lrt_M0_Branch":      "  • M0 (null) vs Branch (alt) - Branch-specific selection\n",
    "log_lrt_BranchSite":     "  • Branch-site_null (null) vs Branch-site (alt) - Branch-site selection\n",
    "log_lrt_total":          "Total of {n} comparisons found.\n\n",
    "log_regenerating":       "Regenerating summary files...\n",
    "log_update_done":        "\n[OK] UPDATE COMPLETED SUCCESSFULLY!\n",
    "log_resolving_stops":    "▶ Resolving all stop codons in this exon automatically...\n",
    "log_analysis_start":     "▶ STARTING BATCH ANALYSIS\n",
    "log_analysis_stopped":   "■ ANALYSIS STOPPED\n",
    "log_analysis_done":      "[OK] ANALYSIS COMPLETE\n",

    # ── App — file dialogs ───────────────────────────────────────────────
    "dialog_select_results_folder": "Select results folder to update summary",

    # ── App — language switch ────────────────────────────────────────────
    "lang_switch_message":  "Restart EasyPAML to apply the language change.",
    "lang_switch_confirm":  "Restart now?",
    "lang_switch_err":      "Could not restart automatically:\n{error}\n\nPlease reopen manually.",

    # ── ResultsViewerWindow — inline messages ────────────────────────────
    "viewer_genes_loaded":      "{n} genes loaded",
    "viewer_gene_not_found":    "Gene not found.",
    "viewer_branch_no_file":    "Branch results file not found.",
    "viewer_branch_read_err":   "Error reading table: {error}",
    "viewer_branch_no_table":   "dN & dS table not found in file.",
    "viewer_branch_invalid":    "Invalid data.",
    "viewer_lrt_parse_err":     "[Error] Error parsing comparison",
    "viewer_sites_subtitle":    "Model: {model}  ·  Analysis: {method}  ·  {omega}",
    "viewer_sites_count":       "{n} sites",

    # ── ResultsViewerWindow — messageboxes ───────────────────────────────
    "msg_warning":              "Warning",
    "msg_success":              "Success",
    "msg_error":                "Error",
    "msg_no_figure":            "No figure to export.",
    "msg_exported_to":          "Exported to:\n{path}",
    "msg_export_err":           "Export error:\n{error}",
    "msg_no_lrt":               "No LRT results found.",
    "msg_excel_exported":       "Excel exported with {n} sheet(s):\n{path}",
    "msg_excel_err":            "Error exporting Excel:\n{error}",
    "msg_csv_exported":         "{n} file(s) exported:\n{files}",
    "msg_no_csv_data":          "No data available to export.",
    "msg_csv_err":              "Error exporting CSV:\n{error}",
    "msg_html_exported":        "HTML report exported:\n{path}",
    "msg_html_err":             "Error exporting HTML:\n{error}",

    # ── ResultsViewerWindow — file dialogs ───────────────────────────────
    "dialog_export_cladogram":  "Export cladogram",
    "dialog_save_csv":          "Save CSVs — choose base name (no extension)",

    # ── Charts (matplotlib) ───────────────────────────────────────────────
    "chart_omega_dist":         "ω Distribution",
    "chart_freq":               "Frequency",
}


# ═══════════════════════════════════════════════════════════════════════════
# PROXY TRANSPARENTE — roteia TEXTS["key"] pelo idioma ativo
# ═══════════════════════════════════════════════════════════════════════════

_AVAILABLE: dict[str, dict] = {
    'pt': TEXTS_PT,
    'en': TEXTS_EN,
}

_current_lang: str = 'en'


def set_language(lang: str) -> None:
    """Ativa o idioma especificado ('pt' ou 'en')."""
    global _current_lang
    if lang in _AVAILABLE:
        _current_lang = lang


def get_language() -> str:
    """Retorna o código do idioma ativo ('pt' ou 'en')."""
    return _current_lang


class _TextProxy:
    """Proxy transparente: TEXTS['key'] sempre lê do idioma ativo.

    Todo código que faz ``from .gui_texts import TEXTS`` continua
    funcionando sem alteração — a troca de idioma é invisível.
    """

    def __getitem__(self, key: str):
        return _AVAILABLE.get(_current_lang, TEXTS_PT)[key]

    def get(self, key: str, default=None):
        return _AVAILABLE.get(_current_lang, TEXTS_PT).get(key, default)

    def __contains__(self, key: str) -> bool:
        return key in _AVAILABLE.get(_current_lang, TEXTS_PT)


TEXTS = _TextProxy()
