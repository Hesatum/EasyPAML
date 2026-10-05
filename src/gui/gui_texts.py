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
    "model_status_default":     "Padrão",
    "model_status_configured":  "Configurado",

    # ── Seção 2 — TreeLabelWindow ────────────────────────────────────────────
    "tree_labeler_sidebar_title": "Instruções",

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
    "section_config":  "CONFIGURAÇÕES",

    "btn_input_folder":  "Pasta de alinhamentos",
    "btn_tree_file":     "Arquivo de árvore",
    "btn_output_folder": "Pasta de Resultados",
    "label_not_selected": "Não selecionado",

    "btn_view_results":          "Ver Resultados",
    "btn_update_results":        "Atualizar Resultados",
    "label_update_results_hint": "Recalcular os arquivos de síntese de uma pasta de resultados",

    "label_omega_initial":    "ω (dN/dS) inicial:",
    "label_remove_gaps":      "Remover colunas com gaps (cleandata = 1)",
    "label_remove_gaps_hint": (
        "cleandata = 1 no codeml (recomendado).\n"
        "Remove as colunas de códons com gap (-), base incerta\n"
        "(N, ?) ou stop codon em qualquer sequência.\n"
        "A tabela de sítios mostra a numeração do SEU alinhamento\n"
        "(coluna 'Pos. alinhamento') e a do codeml, que conta só\n"
        "as colunas que sobraram.\n"
        "Desligue (cleandata = 0) só se souber o que está fazendo."
    ),
    "label_cpus": "CPUs (paralelismo):",

    "label_ignore_stops":      "Ignorar stop codons",
    "label_ignore_stops_hint": (
        "Desligado (padrão): genes com stop codon no meio da\n"
        "sequência NÃO rodam e aparecem como FALHOU, com a\n"
        "sequência e a posição do códon.\n"
        "Ligado: o codeml roda e trata a coluna inteira do stop\n"
        "como dado ausente (com cleandata = 1 ela sai da análise).\n"
        "Um stop no último códon nunca impede a análise."
    ),

    "label_auto_prune":      "Poda automática da árvore",
    "label_auto_prune_hint": (
        "Para cada gene, tira da árvore os táxons que não estão\n"
        "no alinhamento. Sequências do alinhamento que NÃO estão\n"
        "na árvore são EXCLUÍDAS da análise -- o EasyPAML avisa\n"
        "antes de rodar e sugere o nome mais parecido.\n"
        "Desligado: nomes diferentes fazem o codeml falhar."
    ),

    "tab_site_models":  "Modelos de sítio",
    "tab_branch_model": "Modelo de ramos",
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
    "btn_open_output":   "Abrir pasta de resultados",

    "log_header_title":    "LOG DE EXECUÇÃO",
    "log_header_subtitle": "uma mensagem por linha",

    "log_welcome": (
        "EasyPAML -- análise de seleção positiva com PAML/codeml\n"
        "1. Escolha a pasta de alinhamentos (.fasta, .fas, .phy, .phylip)\n"
        "2. Escolha o arquivo de árvore (.nwk, .tree, .tre, .txt)\n"
        "3. Escolha (ou crie) a pasta de resultados\n"
        "4. Ligue os modelos e clique em Iniciar\n\n"
    ),

    "btn_label_branch":     "Marcar ramos (várias marcas)",
    "btn_label_branchsite": "Marcar ramo foreground (branch-site)",

    "model_info_test_type":      "Tipo de Teste:",
    "model_info_params":         "Parâmetros:",
    "model_info_purpose":        "Propósito:",
    "model_info_interpretation": "Interpretação:",
    "model_info_use_case":       "Quando usar:",
    "model_info_references":     "Referências:",

    # ── Seção 4 — ResultsViewerWindow ───────────────────────────────────────
    "viewer_window_title":       "EasyPAML — Painel de Análise",
    "viewer_error_no_tsv":       "Arquivo analysis_summary.tsv não encontrado!",
    "viewer_error_run_analysis": "Execute uma análise para gerar resultados.",

    "viewer_header_title":    "EasyPAML  —  Resultados",
    "viewer_header_subtitle": "Análise de Seleção Positiva  ·  CODEML / PAML",

    "viewer_tab_lrt":                "LRT e p-valores",
    "viewer_tab_omega":              "ω > 1 Global",
    "viewer_tab_sites":              "Sítios sob seleção",
    "viewer_tab_branchsite_classes": "Classes branch-site",
    "viewer_tab_branch":             "Análise de ramos",
    "viewer_tab_export":             "Exportar",

    "stats_total_genes":        "Total de Genes",
    "stats_models_run":         "Modelos Rodados",
    "stats_positive_selection": "Seleção Positiva Global",
    "stats_avg_omega":          "ω Médio",

    "lrt_no_comparisons": "Nenhuma comparação LRT disponível",
    "lrt_label_model":    "Teste:",

    "lrt_branch_warning": (
        "Branch e Branch-site requerem sequências com mais de "
        "200 pb para estimativas confiáveis de ω. Genes mais curtos podem "
        "produzir estimativas instáveis ou não convergir."
    ),

    "lrt_no_data_for_comparison": "Nenhum dado de LRT para esta comparação",

    "lrt_footer_template": (
        "Total: {total} genes  ·  "
        "significativos (q < 0,05, BH): {sig}  ·  "
        "df = {df}  ·  "
        "p em notação científica; q = p corrigido por Benjamini-Hochberg"
    ),

    "pos_sel_tab_title": "Selecao Global — ω > 1 no gene inteiro",
    "pos_sel_tab_criterion": (
        "Critério: ω médio do modelo M2a ou M8 > 1.0  AND  LRT p < 0.05  ·  "
        "Diferente de seleção em sítios específicos (aba Sítios Positivos)"
    ),
    "pos_sel_none_found":      "Nenhum sinal de seleção positiva detectado",
    "pos_sel_criterion_short": "Critério: ω > 1.0  AND  p-valor < 0.05",
    "pos_sel_badge":           "* Positivo",

    "viewer_tab_interpretation": "Interpretação",
    "go_tab_title":              "Candidatos e enriquecimento de GO",
    "go_tab_criterion": (
        "Genes com LRT significativo (M1a×M2a e/ou M7×M8, q < 0.05 corrigido por "
        "Benjamini-Hochberg), ranqueados por significância, cruzados com anotação "
        "funcional GO. Enriquecimento de GO: Fisher exato, candidatos vs. todos os "
        "genes testados, também corrigido por BH."
    ),
    "go_tab_load_button":    "Carregar anotação GO (.tsv)",
    "go_tab_none_loaded":    "Nenhuma anotação carregada",
    "go_tab_none_loaded_sub": "Carregue o TSV de anotação (colunas: gene_id_full, go_biological_process, go_cellular_component, go_molecular_function)",
    "go_tab_no_candidates":  "Nenhum gene com LRT significativo neste resultado",
    "go_tab_enrichment_header": "Termos GO enriquecidos entre os candidatos",
    "go_tab_candidates_header": "Genes candidatos (ranqueados por p-valor)",
    "go_tab_load_error":     "Falha ao carregar/processar a anotação",

    "sites_label_model":    "Modelo:",
    "sites_label_gene":     "Gene:",
    "sites_label_analysis": "Análise:",
    "sites_label_filter":   "Filtrar Pr(w>1) ≥",
    "sites_btn_update":     "Atualizar",

    "sites_file_not_found": "Arquivo não encontrado: {filename}",
    "sites_parse_error":    "Erro ao ler o arquivo:\n{error}",
    "sites_no_sites":       "Nenhum sítio com Pr(ω>1) ≥ {threshold}",

    "sites_table_headers": [
        "Pos. alinhamento",
        "Pos. codeml",
        "AA",
        "Pr(ω>1)",
        "Sig.",
        "ω (média ± EP)",
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
    "branch_no_data_hint":  "Execute o modelo Branch com uma árvore marcada.",
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
    "tree_err_no_tree":       "Nenhuma árvore selecionada.",
    "tree_err_load":          "Erro ao carregar a árvore:\n{error}",

    # ── TreeLabelWindow — diálogos de tag ────────────────────────────────
    "tag_dialog_edit_title":  "Editar Tag",
    "tag_dialog_edit_prompt": "Ramo atual: {tag}\n\nDigite novo número ou 'remover':",
    "tag_dialog_new_title":   "Número da Tag",
    "tag_dialog_new_prompt":  "Digite o número da tag:\n(ex: 1 para #1, 2 para #2)",

    # ── App — CPU label ──────────────────────────────────────────────────
    "label_cpus_detected":    "(detectado: {n} núcleos)",

    # ── App — mensagens de log ───────────────────────────────────────────
    "log_no_tree_selected":   "Selecione um arquivo de árvore (.nwk) primeiro.\n",
    "log_no_output_folder":   "Selecione uma pasta de resultados primeiro.\n",
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
    "log_update_done":        "[OK] Arquivos de síntese atualizados.\n",
    
    "log_analysis_start":     "Iniciando a análise…\n",
    "log_analysis_stopped":   "Análise interrompida.\n",
    "log_analysis_done":      "Análise terminada.\n",

    # ── App — diálogos de arquivo ────────────────────────────────────────
    "dialog_select_results_folder": "Selecione a pasta com resultados para atualizar síntese",

    # ── App — troca de idioma ────────────────────────────────────────────
    "lang_switch_message":  "Reinicie o EasyPAML para aplicar o novo idioma.",
    "lang_switch_confirm":  "Reiniciar agora?",
    "lang_switch_err":      "Não foi possível reiniciar automaticamente:\n{error}\n\nReabra manualmente.",

    # ── ResultsViewerWindow — mensagens inline ───────────────────────────
    "viewer_genes_loaded":      "{n} gene(s) carregado(s)",
    "viewer_gene_not_found":    "Gene não encontrado.",
    "viewer_branch_no_file":    "Arquivo de resultados do modelo Branch não encontrado.",
    "viewer_branch_read_err":   "Erro ao ler tabela: {error}",
    "viewer_branch_no_table":   "Tabela dN & dS não encontrada no arquivo.",
    "viewer_branch_invalid":    "Dados inválidos.",
    "viewer_lrt_parse_err":     "Erro ao ler a comparação",
    "viewer_sites_subtitle":    "Modelo: {model}  ·  Análise: {method}  ·  {omega}",
    "viewer_sites_count":       "{n} sítio(s)",

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

    # ── Novas chaves (0.3.0) ─────────────────────────────────────────────
    "btn_yes": "Sim",
    "btn_no": "Não",
    "btn_about": "Sobre",
    "about_title": "Sobre o EasyPAML",
    "about_text": (
        "EasyPAML {version}\n\n"
        "codeml: {codeml}\n"
        "Versão do codeml: {codeml_version}\n"
        "Python: {python}\n"
        "Sistema: {platform}\n\n"
        "Código e documentação: https://github.com/Hesatum/EasyPAML\n"
        "Métodos (parâmetros do codeml, LRT, correção BH): METODOS.md"
    ),
    "about_codeml_missing": "não encontrado (Linux: sudo apt install paml)",
    "label_codonfreq": "Frequências de códons (CodonFreq):",
    "label_codonfreq_hint": (
        "Modelo de frequências de códons do codeml.\n"
        "2 = F3x4 é o padrão usado na maioria dos estudos com\n"
        "M7/M8 e M1a/M2a. 7 = FMutSel é outro modelo (mais\n"
        "parâmetros) -- se usar, descreva assim nos métodos."
    ),
    "label_ncatg": "Categorias da beta (ncatG, M7/M8):",
    "label_found_alignments": "{n} alinhamento(s) encontrado(s): {names}",
    "label_no_alignments": "Nenhum alinhamento (.fasta, .fas, .phy, .phylip) nesta pasta",
    "label_output_created": "{name} (pasta criada)",
    "hint_select_files": "Escolha a pasta de alinhamentos, o arquivo de árvore e a pasta de resultados para ativar os modelos.",
    "progress_idle": "Nenhuma análise em andamento",
    "progress_template": "Gene {done} de {total}",
    "progress_done": "{done} de {total} genes processados",
    "chk_show_details": "Mostrar detalhes técnicos",
    "preflight_title": "Verificação dos dados",
    "preflight_running": "Verificando os alinhamentos e a árvore…",
    "preflight_heading": "{genes} gene(s) verificado(s): {errors} erro(s) e {warnings} aviso(s)",
    "preflight_explain": (
        "Erros impedem o gene de rodar. Avisos podem mudar o resultado (por exemplo, uma "
        "sequência excluída da análise). Corrija os arquivos e clique em Iniciar de novo, "
        "ou continue mesmo assim."
    ),
    "preflight_tag_error": "ERRO",
    "preflight_tag_warning": "AVISO",
    "preflight_tag_info": "INFO",
    "preflight_general": "(geral)",
    "preflight_btn_fix": "Corrigir e voltar",
    "preflight_btn_continue": "Continuar mesmo assim",
    "preflight_continue_hint": (
        "Continuando, genes com erro aparecem como FALHOU e stop codons são tratados "
        "pelo codeml como dado ausente."
    ),
    "failures_title": "Genes que falharam",
    "failures_text": "{failed} de {total} gene(s) falharam:\n\n{items}\n\nDetalhes em {path}",
    "msg_output_folder_error": "Não foi possível criar a pasta de resultados:\n{error}",
    "dialog_choose_output": "Escolha ou digite o nome de uma pasta nova para os resultados",
    "neutral_window_title": "Modelos nulos e comparações LRT",
    "neutral_header": "Modelos nulos e comparações LRT",
    "neutral_intro": (
        "Quando ativado, o EasyPAML adiciona automaticamente o modelo nulo de cada par "
        "-- não é preciso marcá-lo. Cada comparação usa o teste da razão de verossimilhança "
        "(LRT): 2ΔlnL comparado ao χ² com os graus de liberdade indicados; o p é corrigido "
        "por Benjamini-Hochberg entre os genes (q)."
    ),
    "neutral_footer": "Com a opção ligada, o modelo nulo de cada par selecionado é adicionado automaticamente.",
    "neutral_null_alt": "nulo: {null}  →  alternativo: {alt}",
    "neutral_note": "Nota: {note}",
    "lrt_pairs": [
        ("M1a", "M2a", "#10b981", "M2a vs M1a", "Seleção positiva por sítio",
         "2 graus de liberdade. M1a permite só purificação (ω < 1) e neutralidade (ω = 1); M2a "
         "acrescenta uma classe com ω > 1. Os sítios são identificados por BEB.",
         "O ω₁ = 1 do M1a é imposto pelo próprio codeml (NSsites = 1)."),
        ("M7", "M8", "#ec4899", "M8 vs M7", "Seleção positiva com distribuição beta",
         "2 graus de liberdade. M7: ω segue uma beta entre 0 e 1. M8: beta + uma classe extra "
         "com ω livre. Pode dar significativo só porque há sítios neutros (ω = 1) -- confira o "
         "M8a vs M8.",
         "M7 não usa fix_omega."),
        ("M8a", "M8", "#f472b6", "M8 vs M8a", "Seleção positiva descontando sítios neutros",
         "1 grau de liberdade. M8a é o M8 com a classe extra fixa em ω = 1. Só rejeita o M8a se "
         "houver sítios com ω > 1 (Swanson et al. 2003). p por χ²₁; a mistura 50:50 aparece "
         "no LRT_results.txt só como referência.",
         "M8a usa fix_omega = 1 e omega = 1 -- configurado automaticamente."),
        ("M0", "Branch", "#f59e0b", "Branch vs M0", "Variação de ω entre linhagens",
         "Graus de liberdade = número de grupos de ramos marcados. M0 usa um ω para todos os "
         "ramos; o modelo Branch estima um ω por grupo marcado.",
         "M0 estima ω livremente."),
        ("Branch-site_null", "Branch-site", "#ef4444", "Branch-site vs nulo",
         "Seleção episódica em sítios do ramo foreground (#1)",
         "1 grau de liberdade, p por χ²₁ (recomendação do manual do PAML); a mistura 50:50 "
         "χ²₀/χ²₁ aparece no LRT_results.txt só como referência.",
         "O nulo usa fix_omega = 1 e omega = 1 -- configurado automaticamente."),
    ],
    "model_desc": {
        "M0": "Um único ω para o gene inteiro. Linha de base e nulo do modelo Branch.",
        "M1a": "Purificação (ω < 1) e neutralidade (ω = 1), sem seleção positiva. Nulo do M2a.",
        "M2a": "M1a + classe com ω > 1. M2a vs M1a testa seleção positiva por sítio (BEB).",
        "M7": "ω segue uma beta entre 0 e 1. Nulo do M8.",
        "M8": "Beta + classe com ω livre. Testado contra M7 e contra M8a.",
        "M8a": "M8 com a classe extra fixa em ω = 1. Nulo do M8 que aceita sítios neutros.",
        "Branch": "ω por grupo de ramos marcado. Testado contra M0.",
        "Branch-site": "Seleção positiva em sítios do ramo foreground (#1).",
        "Branch-site_null": "Branch-site com ω₂ = 1 fixo. Nulo do Branch-site.",
    },
    "cfg_btn": "editar",
    "cfg_field_nssites": ("NSsites", "Modelo de sítios: 0 = M0, 1 = M1a, 2 = M2a, 7 = M7, 8 = M8"),
    "cfg_field_model": ("model", "0 = um ω para todos os ramos; 2 = ω por grupo de ramos"),
    "cfg_field_fix_omega": ("fix_omega", "0 = ω estimado; 1 = ω fixo no valor abaixo"),
    "cfg_field_omega": ("omega", "ω inicial (ou fixo, se fix_omega = 1)"),
    "cfg_field_codonfreq": ("CodonFreq", "Modelo de frequências de códons"),
    "cfg_field_ncatg": ("ncatG", "Categorias da distribuição beta (M7/M8/M8a)"),
    "cfg_field_kappa": ("kappa", "κ (ts/tv) inicial"),
    "cfg_preview": "Parâmetros que irão para o .ctl deste modelo:",
    "cfg_saved": "[OK] Parâmetros de {model} salvos.",
    "lang_restart_title": "Idioma",
    "viewer_tab_summary": "Resumo",
    "viewer_btn_open_output": "Abrir pasta de resultados",
    "stats_sig_genes": "Genes com seleção positiva (q < 0,05)",
    "stats_failed": "Genes que falharam",
    "summary_title": "Uma frase por gene",
    "summary_explain": (
        "Para cada teste: p do LRT, q (p corrigido por Benjamini-Hochberg entre os genes), "
        "ω e proporção (p₁) da classe de sítios que pode ter ω > 1, e quantos sítios têm "
        "Pr(ω>1) ≥ 0,95 no BEB. O ω médio do gene NÃO é critério de seleção positiva: ele "
        "fica abaixo de 1 mesmo quando poucos sítios estão sob seleção forte."
    ),
    "summary_sig": "{test}: significativo (p = {p}, q = {q}) — {sites} sítio(s) com Pr(ω>1) ≥ 0,95",
    "summary_nonsig": "{test}: não significativo (p = {p}, q = {q})",
    "summary_posclass": "classe positiva: ω = {w}, p₁ = {p1}",
    "summary_failed": "FALHOU: {reason}",
    "summary_m8a_caveat": (
        "Atenção: M7 vs M8 é significativo mas M8a vs M8 não -- o sinal pode vir de sítios "
        "neutros (ω = 1), não de seleção positiva."
    ),
    "summary_no_tests": "Nenhum teste de seleção positiva (M1a/M2a, M7/M8, M8a/M8) neste resultado.",
    "sites_legend": "* Pr(ω>1) ≥ 0,95   ** Pr(ω>1) ≥ 0,99   (Pr = probabilidade posterior de o sítio estar na classe com ω > 1)",
    "sites_position_note": (
        "'Pos. alinhamento' é a posição do códon no SEU arquivo. 'Pos. codeml' é a do arquivo "
        "bruto do codeml, que conta só as colunas que sobraram depois da limpeza (cleandata)."
    ),
    "sites_unmapped_note": (
        "Resultado sem mapa de numeração (gerado por versão antiga ou sem verificação): "
        "as duas colunas mostram a numeração do codeml."
    ),
    "sites_btn_copy": "Copiar sítios (TSV)",
    "sites_btn_export": "Exportar TSV…",
    "sites_copied": "{n} sítio(s) copiado(s) para a área de transferência.",
    "lrt_headers": ["Gene", "lnL nulo", "lnL alternativo", "2Δℓ", "p", "q (BH)", "ω classe + (p₁)", "Sig."],
    "lrt_sig_yes": "sim",
    "lrt_sig_no": "não",
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

    "label_omega_initial":    "Initial ω (dN/dS):",
    "label_remove_gaps":      "Remove gap columns (cleandata = 1)",
    "label_remove_gaps_hint": (
        "cleandata = 1 in codeml (recommended).\n"
        "Removes codon columns with a gap (-), an ambiguous base\n"
        "(N, ?) or a stop codon in any sequence.\n"
        "The sites table shows YOUR alignment numbering\n"
        "('Aln. pos.' column) and codeml's, which only counts\n"
        "the remaining columns.\n"
        "Turn off (cleandata = 0) only if you know why."
    ),
    "label_cpus": "CPUs (parallelism):",

    "label_ignore_stops":      "Ignore stop codons",
    "label_ignore_stops_hint": (
        "Off (default): genes with a stop codon inside the\n"
        "sequence do NOT run and are shown as FAILED, with the\n"
        "sequence and codon position.\n"
        "On: codeml runs and treats the whole stop column as\n"
        "missing data (with cleandata = 1 it is removed).\n"
        "A stop in the last codon never blocks the analysis."
    ),

    "label_auto_prune":      "Automatic tree pruning",
    "label_auto_prune_hint": (
        "For each gene, removes from the tree the taxa that are\n"
        "not in the alignment. Alignment sequences that are NOT\n"
        "in the tree are EXCLUDED from the analysis -- EasyPAML\n"
        "warns before running and suggests the closest name.\n"
        "Off: mismatched names make codeml fail."
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
    "btn_open_output":   "Open results folder",

    "log_header_title":    "EXECUTION LOG",
    "log_header_subtitle": "one message per line",

    "log_welcome": (
        "EasyPAML -- positive selection analysis with PAML/codeml\n"
        "1. Choose the alignments folder (.fasta, .fas, .phy, .phylip)\n"
        "2. Choose the tree file (.nwk, .tree, .tre, .txt)\n"
        "3. Choose (or create) the results folder\n"
        "4. Switch models on and click Run\n\n"
    ),

    "btn_label_branch":     "Label branches (multiple tags)",
    "btn_label_branchsite": "Label foreground branch (branch-site)",

    "model_info_test_type":      "Test Type:",
    "model_info_params":         "Parameters:",
    "model_info_purpose":        "Purpose:",
    "model_info_interpretation": "Interpretation:",
    "model_info_use_case":       "When to use:",
    "model_info_references":     "References:",

    # ── Section 4 — ResultsViewerWindow ─────────────────────────────────────
    "viewer_window_title":       "EasyPAML — Analysis Panel",
    "viewer_error_no_tsv":       "File analysis_summary.tsv not found!",
    "viewer_error_run_analysis": "Run an analysis to generate results.",

    "viewer_header_title":    "EasyPAML  —  Results",
    "viewer_header_subtitle": "Positive Selection Analysis  ·  CODEML / PAML",

    "viewer_tab_lrt":                "LRT and p-values",
    "viewer_tab_omega":              "ω > 1 Global",
    "viewer_tab_sites":              "Positive Sites",
    "viewer_tab_branchsite_classes": "Branch-site Classes",
    "viewer_tab_branch":             "Branch Analysis",
    "viewer_tab_export":             "Export",

    "stats_total_genes":        "Total Genes",
    "stats_models_run":         "Models Run",
    "stats_positive_selection": "Global Positive Selection",
    "stats_avg_omega":          "Average ω",

    "lrt_no_comparisons": "No LRT comparison available",
    "lrt_label_model":    "Test:",

    "lrt_branch_warning": (
        "Branch and Branch-site require sequences longer than "
        "200 bp for reliable ω estimates. Shorter genes may "
        "produce unstable estimates or fail to converge."
    ),

    "lrt_no_data_for_comparison": "No LRT data for this comparison",

    "lrt_footer_template": (
        "Total: {total} genes  ·  "
        "significant (q < 0.05, BH): {sig}  ·  "
        "df = {df}  ·  "
        "p in scientific notation; q = Benjamini-Hochberg corrected p"
    ),

    "pos_sel_tab_title": "Global Selection — ω > 1 across the whole gene",
    "pos_sel_tab_criterion": (
        "Criterion: mean ω from M2a or M8 > 1.0  AND  LRT p < 0.05  ·  "
        "Different from site-specific selection (Positive Sites tab)"
    ),
    "pos_sel_none_found":      "No positive selection signal detected",
    "pos_sel_criterion_short": "Criterion: ω > 1.0  AND  p-value < 0.05",
    "pos_sel_badge":           "* Positive",

    "viewer_tab_interpretation": "Interpretation",
    "go_tab_title":              "Candidates and GO enrichment",
    "go_tab_criterion": (
        "Genes with significant LRT (M1a×M2a and/or M7×M8, Benjamini-Hochberg "
        "q < 0.05), ranked by significance, cross-referenced with GO functional "
        "annotation. GO enrichment: Fisher exact test, candidates vs. all tested "
        "genes, also BH-corrected."
    ),
    "go_tab_load_button":    "Load GO annotation (.tsv)",
    "go_tab_none_loaded":    "No annotation loaded",
    "go_tab_none_loaded_sub": "Load the annotation TSV (columns: gene_id_full, go_biological_process, go_cellular_component, go_molecular_function)",
    "go_tab_no_candidates":  "No gene with significant LRT in this result",
    "go_tab_enrichment_header": "GO terms enriched among candidates",
    "go_tab_candidates_header": "Candidate genes (ranked by p-value)",
    "go_tab_load_error":     "Failed to load/process the annotation",

    "sites_label_model":    "Model:",
    "sites_label_gene":     "Gene:",
    "sites_label_analysis": "Analysis:",
    "sites_label_filter":   "Filter Pr(w>1) ≥",
    "sites_btn_update":     "Update",

    "sites_file_not_found": "File not found: {filename}",
    "sites_parse_error":    "Error reading the file:\n{error}",
    "sites_no_sites":       "No site with Pr(ω>1) ≥ {threshold}",

    "sites_table_headers": [
        "Aln. pos.",
        "codeml pos.",
        "AA",
        "Pr(ω>1)",
        "Sig.",
        "ω (mean ± SE)",
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
    "tree_err_no_tree":       "No tree selected.",
    "tree_err_load":          "Failed to load the tree:\n{error}",

    # ── TreeLabelWindow — tag dialogs ────────────────────────────────────
    "tag_dialog_edit_title":  "Edit Tag",
    "tag_dialog_edit_prompt": "Current branch: {tag}\n\nEnter new number or 'remove':",
    "tag_dialog_new_title":   "Tag Number",
    "tag_dialog_new_prompt":  "Enter tag number:\n(e.g. 1 for #1, 2 for #2)",

    # ── App — CPU label ──────────────────────────────────────────────────
    "label_cpus_detected":    "(detected: {n} cores)",

    # ── App — log messages ───────────────────────────────────────────────
    "log_no_tree_selected":   "Please choose a tree file (.nwk) first.\n",
    "log_no_output_folder":   "Please choose a results folder first.\n",
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
    "log_update_done":        "[OK] Summary files updated.\n",
    
    "log_analysis_start":     "Starting the analysis…\n",
    "log_analysis_stopped":   "Analysis stopped.\n",
    "log_analysis_done":      "Analysis finished.\n",

    # ── App — file dialogs ───────────────────────────────────────────────
    "dialog_select_results_folder": "Select results folder to update summary",

    # ── App — language switch ────────────────────────────────────────────
    "lang_switch_message":  "Restart EasyPAML to apply the language change.",
    "lang_switch_confirm":  "Restart now?",
    "lang_switch_err":      "Could not restart automatically:\n{error}\n\nPlease reopen manually.",

    # ── ResultsViewerWindow — inline messages ────────────────────────────
    "viewer_genes_loaded":      "{n} gene(s) loaded",
    "viewer_gene_not_found":    "Gene not found.",
    "viewer_branch_no_file":    "Branch results file not found.",
    "viewer_branch_read_err":   "Error reading table: {error}",
    "viewer_branch_no_table":   "dN & dS table not found in file.",
    "viewer_branch_invalid":    "Invalid data.",
    "viewer_lrt_parse_err":     "[Error] Error parsing comparison",
    "viewer_sites_subtitle":    "Model: {model}  ·  Analysis: {method}  ·  {omega}",
    "viewer_sites_count":       "{n} site(s)",

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

    # ── New keys (0.3.0) ─────────────────────────────────────────────────
    "btn_yes": "Yes",
    "btn_no": "No",
    "btn_about": "About",
    "about_title": "About EasyPAML",
    "about_text": (
        "EasyPAML {version}\n\n"
        "codeml: {codeml}\n"
        "codeml version: {codeml_version}\n"
        "Python: {python}\n"
        "System: {platform}\n\n"
        "Code and documentation: https://github.com/Hesatum/EasyPAML\n"
        "Methods (codeml parameters, LRT, BH correction): METODOS.md"
    ),
    "about_codeml_missing": "not found (Linux: sudo apt install paml)",
    "label_codonfreq": "Codon frequencies (CodonFreq):",
    "label_codonfreq_hint": (
        "codeml codon frequency model.\n"
        "2 = F3x4 is the default used in most M7/M8 and\n"
        "M1a/M2a studies. 7 = FMutSel is a different model\n"
        "(more parameters) -- if you use it, say so in Methods."
    ),
    "label_ncatg": "Beta categories (ncatG, M7/M8):",
    "label_found_alignments": "{n} alignment(s) found: {names}",
    "label_no_alignments": "No alignment (.fasta, .fas, .phy, .phylip) in this folder",
    "label_output_created": "{name} (folder created)",
    "hint_select_files": "Choose the alignments folder, the tree file and the results folder to enable the models.",
    "progress_idle": "No analysis running",
    "progress_template": "Gene {done} of {total}",
    "progress_done": "{done} of {total} genes processed",
    "chk_show_details": "Show technical details",
    "preflight_title": "Data check",
    "preflight_running": "Checking alignments and tree…",
    "preflight_heading": "{genes} gene(s) checked: {errors} error(s) and {warnings} warning(s)",
    "preflight_explain": (
        "Errors prevent a gene from running. Warnings may change the result (for example, a "
        "sequence excluded from the analysis). Fix the files and click Run again, or continue anyway."
    ),
    "preflight_tag_error": "ERROR",
    "preflight_tag_warning": "WARNING",
    "preflight_tag_info": "INFO",
    "preflight_general": "(general)",
    "preflight_btn_fix": "Fix and go back",
    "preflight_btn_continue": "Continue anyway",
    "preflight_continue_hint": (
        "If you continue, genes with errors are shown as FAILED and stop codons are "
        "treated by codeml as missing data."
    ),
    "failures_title": "Genes that failed",
    "failures_text": "{failed} of {total} gene(s) failed:\n\n{items}\n\nDetails in {path}",
    "msg_output_folder_error": "Could not create the results folder:\n{error}",
    "dialog_choose_output": "Choose or type the name of a new folder for the results",
    "neutral_window_title": "Null models and LRT comparisons",
    "neutral_header": "Null models and LRT comparisons",
    "neutral_intro": (
        "When enabled, EasyPAML automatically adds the null model of each pair -- you don't "
        "need to select it. Each comparison uses the likelihood ratio test (LRT): 2ΔlnL "
        "compared to χ² with the degrees of freedom shown; p is Benjamini-Hochberg corrected "
        "across genes (q)."
    ),
    "neutral_footer": "With this option on, the null model of each selected pair is added automatically.",
    "neutral_null_alt": "null: {null}  →  alternative: {alt}",
    "neutral_note": "Note: {note}",
    "lrt_pairs": [
        ("M1a", "M2a", "#10b981", "M2a vs M1a", "Site-wise positive selection",
         "2 degrees of freedom. M1a allows only purifying (ω < 1) and neutral (ω = 1) sites; M2a "
         "adds a class with ω > 1. Sites are identified by BEB.",
         "M1a's ω₁ = 1 is imposed by codeml itself (NSsites = 1)."),
        ("M7", "M8", "#ec4899", "M8 vs M7", "Positive selection with a beta distribution",
         "2 degrees of freedom. M7: ω follows a beta between 0 and 1. M8: beta + one extra class "
         "with free ω. It can be significant just because some sites are neutral (ω = 1) -- check "
         "M8a vs M8.",
         "M7 does not use fix_omega."),
        ("M8a", "M8", "#f472b6", "M8 vs M8a", "Positive selection beyond neutral sites",
         "1 degree of freedom. M8a is M8 with the extra class fixed at ω = 1. It only rejects M8a "
         "if there are sites with ω > 1 (Swanson et al. 2003). p from χ²₁; the 50:50 mixture is "
         "reported in LRT_results.txt as a reference only.",
         "M8a uses fix_omega = 1 and omega = 1 -- set automatically."),
        ("M0", "Branch", "#f59e0b", "Branch vs M0", "ω variation among lineages",
         "Degrees of freedom = number of labelled branch groups. M0 uses one ω for all branches; "
         "the Branch model estimates one ω per labelled group.",
         "M0 estimates ω freely."),
        ("Branch-site_null", "Branch-site", "#ef4444", "Branch-site vs null",
         "Episodic selection at sites of the foreground branch (#1)",
         "1 degree of freedom, p from χ²₁ (PAML manual recommendation); the 50:50 χ²₀/χ²₁ "
         "mixture is reported in LRT_results.txt as a reference only.",
         "The null uses fix_omega = 1 and omega = 1 -- set automatically."),
    ],
    "model_desc": {
        "M0": "A single ω for the whole gene. Baseline and null of the Branch model.",
        "M1a": "Purifying (ω < 1) and neutral (ω = 1) sites, no positive selection. Null of M2a.",
        "M2a": "M1a + a class with ω > 1. M2a vs M1a tests site-wise positive selection (BEB).",
        "M7": "ω follows a beta between 0 and 1. Null of M8.",
        "M8": "Beta + a class with free ω. Tested against M7 and against M8a.",
        "M8a": "M8 with the extra class fixed at ω = 1. Null of M8 that allows neutral sites.",
        "Branch": "One ω per labelled branch group. Tested against M0.",
        "Branch-site": "Positive selection at sites of the foreground branch (#1).",
        "Branch-site_null": "Branch-site with ω₂ = 1 fixed. Null of Branch-site.",
    },
    "cfg_btn": "edit",
    "cfg_field_nssites": ("NSsites", "Site model: 0 = M0, 1 = M1a, 2 = M2a, 7 = M7, 8 = M8"),
    "cfg_field_model": ("model", "0 = one ω for all branches; 2 = ω per branch group"),
    "cfg_field_fix_omega": ("fix_omega", "0 = ω estimated; 1 = ω fixed at the value below"),
    "cfg_field_omega": ("omega", "initial ω (fixed value if fix_omega = 1)"),
    "cfg_field_codonfreq": ("CodonFreq", "Codon frequency model"),
    "cfg_field_ncatg": ("ncatG", "Categories of the beta distribution (M7/M8/M8a)"),
    "cfg_field_kappa": ("kappa", "initial κ (ts/tv)"),
    "cfg_preview": "Parameters that will go into this model's .ctl:",
    "cfg_saved": "[OK] {model} parameters saved.",
    "lang_restart_title": "Language",
    "viewer_tab_summary": "Summary",
    "viewer_btn_open_output": "Open results folder",
    "stats_sig_genes": "Genes with positive selection (q < 0.05)",
    "stats_failed": "Genes that failed",
    "summary_title": "One sentence per gene",
    "summary_explain": (
        "For each test: LRT p, q (p corrected by Benjamini-Hochberg across genes), ω and "
        "proportion (p₁) of the site class that can have ω > 1, and how many sites have "
        "Pr(ω>1) ≥ 0.95 in BEB. The gene's mean ω is NOT a positive selection criterion: it "
        "stays below 1 even when a few sites are under strong selection."
    ),
    "summary_sig": "{test}: significant (p = {p}, q = {q}) — {sites} site(s) with Pr(ω>1) ≥ 0.95",
    "summary_nonsig": "{test}: not significant (p = {p}, q = {q})",
    "summary_posclass": "positive class: ω = {w}, p₁ = {p1}",
    "summary_failed": "FAILED: {reason}",
    "summary_m8a_caveat": (
        "Caution: M7 vs M8 is significant but M8a vs M8 is not -- the signal may come from "
        "neutral sites (ω = 1), not from positive selection."
    ),
    "summary_no_tests": "No positive selection test (M1a/M2a, M7/M8, M8a/M8) in this result.",
    "sites_legend": "* Pr(ω>1) ≥ 0.95   ** Pr(ω>1) ≥ 0.99   (Pr = posterior probability that the site is in the ω > 1 class)",
    "sites_position_note": (
        "'Aln. pos.' is the codon position in YOUR file. 'codeml pos.' is the one in codeml's raw "
        "output, which only counts the columns left after cleaning (cleandata)."
    ),
    "sites_unmapped_note": (
        "Result without a numbering map (from an older version or unverified): both columns "
        "show codeml's numbering."
    ),
    "sites_btn_copy": "Copy sites (TSV)",
    "sites_btn_export": "Export TSV…",
    "sites_copied": "{n} site(s) copied to the clipboard.",
    "lrt_headers": ["Gene", "lnL null", "lnL alternative", "2Δℓ", "p", "q (BH)", "positive-class ω (p₁)", "Sig."],
    "lrt_sig_yes": "yes",
    "lrt_sig_no": "no",
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


def tr(pt: str, en: str) -> str:
    """Mensagem curta bilíngue (para textos montados em código)."""
    return pt if _current_lang == 'pt' else en


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
