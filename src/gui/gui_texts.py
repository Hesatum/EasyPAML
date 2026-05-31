"""
gui_texts.py — Dicionário Central de Textos da Interface EasyPAML
=================================================================

Todos os textos visíveis ao usuário estão centralizados aqui.
Para editar qualquer rótulo, botão ou mensagem da interface, basta
modificar o valor correspondente neste arquivo — não é necessário
tocar em main_gui.py ou results_viewer.py.

ORGANIZAÇÃO
-----------
  Seção 1 — main_gui.py › ModelConfigWindow   (janela de configuração do modelo)
  Seção 2 — main_gui.py › TreeLabelWindow      (janela de etiquetagem da árvore)
  Seção 3 — main_gui.py › App                  (janela principal, sidebar, log)
  Seção 4 — results_viewer.py › ResultsViewerWindow (painel de resultados)

TEMPLATES (strings com {placeholders})
---------------------------------------
  Algumas strings contêm variáveis dinâmicas entre chaves.
  Use .format() para substituí-las. Exemplo:
    TEXTS["sites_file_not_found"].format(filename="gene_M8.txt")
    TEXTS["lrt_footer_template"].format(total=30, sig=5, df=2)
    TEXTS["status_stops_template"].format(n=3)
"""

TEXTS: dict[str, object] = {

    # ═══════════════════════════════════════════════════════════════════════
    # SEÇÃO 1 — main_gui.py › ModelConfigWindow
    # Janela de configuração dos parâmetros .ctl de cada modelo CODEML.
    # Aberta clicando em "cfg" no card do modelo.
    # ═══════════════════════════════════════════════════════════════════════

    # Cabeçalho da janela — exibe o código do modelo (ex: "Modelo: M8")
    "model_config_header": "Modelo: {model_code}",

    # Botões do rodapé da janela de configuração
    "model_config_btn_cancel": "Cancelar",
    "model_config_btn_save":   "Salvar",

    # Status exibido no card do modelo na janela principal
    # "padrao"      → antes de qualquer configuração manual
    # "Configurado" → após salvar parâmetros personalizados
    "model_status_default":    "Default",
    "model_status_configured": "Configurado",


    # ═══════════════════════════════════════════════════════════════════════
    # SEÇÃO 2 — main_gui.py › TreeLabelWindow
    # Janela interativa de etiquetagem de ramos da árvore filogenética.
    # Aberta pelos botões "Marcar Ramos" e "Marcar Branch-site".
    # ═══════════════════════════════════════════════════════════════════════

    # Título da barra lateral esquerda da janela de etiquetagem
    "tree_labeler_sidebar_title": "Instrucoes",

    # Instruções exibidas quando modo = 'branchsite'
    # (apenas a tag #1/foreground é permitida)
    "tree_labeler_instructions_branchsite": (
        "• Clique nos círculos\n  para marcar/desmarcar\n  foreground (#1).\n\n"
        "• O modo branch-site só permite uma tag #1\n.\n\n"
        "• Cor: Vermelho"
    ),

    # Instruções exibidas quando modo = 'branch'
    # (múltiplas tags numéricas, cada uma com cor única)
    "tree_labeler_instructions_branch": (
        "• Clique nos círculos\n  para atribuir tags e testar o ômega dos ramos.\n\n"
        "• Digite número\n  da tag (1, 2, 3...).\n\n"
        "• Cada número\n  recebe cor única."
    ),

    # Legenda de tags ativas na sidebar
    "tree_labeler_legend_title": "Tags Ativas",
    # Exibido quando nenhum ramo foi etiquetado ainda
    "tree_labeler_no_tags":      "Nenhuma tag ativa",

    # Botões SALVAR / CANCELAR da TreeLabelWindow
    "tree_labeler_btn_save":   "Salvar",
    "tree_labeler_btn_cancel": "Cancelar",


    # ═══════════════════════════════════════════════════════════════════════
    # SEÇÃO 3 — main_gui.py › App (janela principal)
    # ═══════════════════════════════════════════════════════════════════════

    # ── Cabeçalho da sidebar ─────────────────────────────────────────────
    # Logo e subtítulo no topo da barra lateral esquerda
    "app_sidebar_title":    "EasyPAML",
    "app_sidebar_subtitle": "Seleção Positiva",

    # ── Títulos das seções colapsáveis da sidebar ────────────────────────
    "section_files":   "ARQUIVOS",
    "section_results": "RESULTADOS",
    "section_config":  "CONFIGURACOES",

    # ── Seção ARQUIVOS — botões e labels de seleção ──────────────────────
    # Texto dos três botões de seleção de caminho
    "btn_input_folder":  "Pasta de alinhamentos",
    "btn_tree_file":     "Tree file",
    "btn_output_folder": "Pasta de Resultados",
    # Label exibido abaixo de cada botão antes de qualquer seleção
    "label_not_selected": "Nao selecionado",

    # ── Seção RESULTADOS — botões ────────────────────────────────────────
    "btn_view_results":          "Ver Resultados",
    "btn_update_results":        "Atualizar Resultados",
    # Dica exibida abaixo do botão de atualização
    "label_update_results_hint": "Atualizar arquivos de analise",

    # ── Seção CONFIGURACOES — labels dos campos ──────────────────────────
    "label_omega_initial":  "dN/dS Inicial (w):",
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
    "label_cpus":           "CPUs (paralelismo):",

    # Toggle de modo heurístico (fixa κ do M0 nos modelos seguintes)
    "label_heuristic_mode": "Modo Heuristico (acelerar analise)",
    # Dica textual abaixo do toggle
    "label_heuristic_hint": (
        "Roda M0 como pré-passo e usa seus resultados\n"
        "para acelerar os demais modelos de duas formas:\n"
        "  1. Branch lengths do M0 como ponto de partida\n"
        "     (maior ganho — evita busca do zero).\n"
        "  2. Kappa (ts/tv) do M0 fixado nos modelos\n"
        "     seguintes (ganho menor, ~5-15%).\n"
        "Se M0 já estiver selecionado, o warm-start de\n"
        "branch lengths ocorre automaticamente; este\n"
        "toggle adiciona apenas o kappa fixo."
    ),

    # Toggle de modo WGS (ativa ndata no .ctl)
    "label_wgs_mode": "Modo lote unico (ndata)",
    "label_wgs_hint": (
        "Combina TODOS os genes em um unico arquivo e\n"
        "executa o CODEML uma vez com ndata=N.\n"
        "Exige que todos os genes tenham EXATAMENTE\n"
        "as mesmas especies (sem poda automatica).\n"
        "Use apenas se os genes sao homogeneos.\n"
        "Para analise padrao por gene, deixe desativado."
    ),

    # Toggle de ignorar stop codons (auto_continue_stop_codons)
    # Desativado por padrão — CODEML irá parar e reportar stop codons normalmente.
    # Quando ativado, envia Enter automaticamente sem interromper a análise.
    # No modo paralelo (CPUs > 1) este comportamento é sempre forçado.
    "label_ignore_stops":      "Ignorar Stop Codons",
    "label_ignore_stops_hint": (
        "Envia Enter automaticamente ao encontrar stop codons.\n"
        "No modo paralelo (CPUs > 1) é sempre ativado."
    ),

    # Poda automática de árvore
    "label_auto_prune":      "Poda Automática de Árvore",
    "label_auto_prune_hint": (
        "Poda a árvore para incluir apenas os\n"
        "taxons presentes em cada gene automaticamente."
    ),

    # ── Abas dos modelos (CTkTabview central) ────────────────────────────
    "tab_site_models":  "Site Models",
    "tab_branch_model": "Branch Model",
    "tab_branchsite":   "Branch-Site",

    # ── Barra de status (linha superior da área de controle) ─────────────
    # Estado inicial — análise não iniciada
    "status_ready": "● Pronto",
    # Estados dinâmicos atualizados durante a análise
    "status_running": "● Executando",
    "status_paused":  "● Pausada",
    "status_stopped": "● Parada",
    # Template do contador de codons STOP detectados
    # Uso: TEXTS["status_stops_template"].format(n=3)
    "status_stops_template": "■ Stops: {n}",
    # Label do toggle de modelos nulos automáticos
    "label_neutral_models": "Modelos nulos automáticos",

    # ── Botões de controle de análise ────────────────────────────────────
    "btn_run":           "Iniciar",
    "btn_pause":         "Pausar",
    # Texto do botão de pause quando a análise está pausada (toggle)
    "btn_resume":        "Retomar",
    "btn_stop":          "Parar",
    "btn_resolve_stops": "Ignorar Stop codons",

    # ── Log de execução ──────────────────────────────────────────────────
    # Rótulo e subtítulo do cabeçalho acima da caixa de texto
    "log_header_title":    "LOG DE EXECUCAO",
    "log_header_subtitle": "output em tempo real",

    # Mensagem de boas-vindas inserida ao iniciar o app
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

    # ── Botões de etiquetagem da árvore (dentro das abas de modelos) ─────
    "btn_label_branch":     "Marcar Branch (Multiplos Ramos)",
    "btn_label_branchsite": "Marcar Branch-site",

    # ── Janela de detalhes de cada modelo (_show_model_info) ─────────────
    # Rótulos das seções de informação (abertas pelo botão "?" no card)
    "model_info_test_type":      "Tipo de Teste:",
    "model_info_params":         "Parametros:",
    "model_info_purpose":        "Proposito:",
    "model_info_interpretation": "Interpretacao:",
    "model_info_use_case":       "Quando usar:",
    "model_info_references":     "Referencias:",


    # ═══════════════════════════════════════════════════════════════════════
    # SEÇÃO 4 — results_viewer.py › ResultsViewerWindow
    # Painel de análise de resultados CODEML, aberto por "Ver Resultados".
    # ═══════════════════════════════════════════════════════════════════════

    # ── Janela principal ─────────────────────────────────────────────────
    # Título exibido na barra de título do sistema operacional
    "viewer_window_title": "EasyPAML - Painel de Analise",
    # Mensagem de erro quando analysis_summary.tsv não é encontrado
    "viewer_error_no_tsv":       "Arquivo analysis_summary.tsv não encontrado!",
    "viewer_error_run_analysis": "Execute uma análise para gerar resultados.",

    # ── Cabeçalho (header strip horizontal) ──────────────────────────────
    "viewer_header_title":    "EasyPAML  —  Resultados",
    "viewer_header_subtitle": "Análise de Seleção Positiva  ·  CODEML / PAML",

    # ── Abas do painel de resultados ─────────────────────────────────────
    "viewer_tab_lrt":                "LRT & P-values",
    "viewer_tab_omega":              "ω > 1 Global",
    "viewer_tab_sites":              "Sítios Positivos",
    "viewer_tab_branchsite_classes": "Branch-site Classes",
    "viewer_tab_branch":             "Branch Analysis",
    "viewer_tab_export":             "Exportar",

    # ── Painel de estatísticas (quatro cards superiores) ─────────────────
    "stats_total_genes":        "Total de Genes",
    "stats_models_run":         "Modelos Rodados",
    "stats_positive_selection": "Seleção Positiva Global",
    "stats_avg_omega":          "ω Médio",

    # ── Aba LRT & P-values ────────────────────────────────────────────────
    # Exibido quando nenhuma coluna lrt_* está no TSV
    "lrt_no_comparisons": "[!] Sem comparacoes LRT disponiveis",
    # Label do combo de seleção de comparação
    "lrt_label_model":    "Modelo testado:",

    # Aviso exibido para comparações Branch/Branch-site (requerem > 200 pb)
    "lrt_branch_warning": (
        "Branch e Branch-site requerem sequencias com mais de "
        "200 pb para estimativas confiaveis de ω.  Genes mais curtos podem "
        "produzir estimativas instaveis ou nao convergir."
    ),

    # Exibido quando nenhum gene tem dados LRT para a comparação selecionada
    "lrt_no_data_for_comparison": "[!] Nenhum dado LRT disponivel para esta comparacao",

    # Rodapé da tabela LRT com contagem de genes significantes
    # Uso: TEXTS["lrt_footer_template"].format(total=30, sig=5, df=2)
    "lrt_footer_template": (
        "Total: {total} genes  ·  "
        "Significantes p < 0.05: {sig}  ·  "
        "df = {df}  ·  "
        "* = significante  ·  * ω>1 = selecao positiva"
    ),

    # ── Aba Global ω > 1 (Seleção Positiva Global) ───────────────────────
    # Título e critério exibidos no banner informativo do topo
    "pos_sel_tab_title": "Selecao Global — ω > 1 no gene inteiro",
    "pos_sel_tab_criterion": (
        "Critério: ω médio do modelo M2a ou M8 > 1.0  AND  LRT p < 0.05  ·  "
        "Diferente de seleção em sítios específicos (aba Sítios Positivos)"
    ),
    # Exibido quando nenhum gene passa no critério de seleção positiva
    "pos_sel_none_found":      "Nenhum sinal de seleção positiva detectado",
    "pos_sel_criterion_short": "Critério: ω > 1.0  AND  p-valor < 0.05",
    # Texto do badge verde em cada gene com seleção detectada
    "pos_sel_badge": "* Positivo",

    # ── Aba Positive Sites ────────────────────────────────────────────────
    # Labels dos combos de seleção de modelo/gene/análise
    "sites_label_model":    "Modelo:",
    "sites_label_gene":     "Gene:",
    "sites_label_analysis": "Análise:",
    # Label do filtro de probabilidade posterior
    "sites_label_filter":   "Filtrar Pr(w>1) ≥",
    # Botão de atualização da tabela
    "sites_btn_update":     "Atualizar",

    # Mensagens de erro/aviso — use .format() para preencher os placeholders
    # Uso: TEXTS["sites_file_not_found"].format(filename="gene_M8_results.txt")
    "sites_file_not_found": "[X] Arquivo nao encontrado: {filename}",
    # Uso: TEXTS["sites_parse_error"].format(error=str(e))
    "sites_parse_error":    "[Erro] Erro ao parsear arquivo:\n{error}",
    # Uso: TEXTS["sites_no_sites"].format(threshold=0.95)
    "sites_no_sites":       "Nenhum sitio acima do limiar Pr(w>1) >= {threshold}",

    # Cabeçalhos das colunas da tabela de sítios (ordem fixa)
    # As larguras em pixels ficam no código — edite apenas os rótulos aqui
    "sites_table_headers": [
        "Posição",    # largura 70 px
        "AA",         # largura 55 px
        "Pr(ω>1)",    # largura 100 px
        "Sig",        # largura 55 px
        "ω (média)",  # largura 105 px
        "ω − SE",     # largura 105 px
        "ω + SE",     # largura 105 px
    ],

    # ── Aba Branch-site Classes ───────────────────────────────────────────
    # Label do combo de seleção de gene
    "branchsite_classes_gene_label": "Gene:",
    # Exibido quando o gene digitado não é encontrado no DataFrame
    "branchsite_classes_not_found":  "Gene não encontrado",
    # Título do bloco de classes — {gene} = nome do gene selecionado
    "branchsite_classes_header":     "Classes do Modelo Branch-site — {gene}",

    # Cabeçalhos das colunas da tabela de classes (ordem fixa)
    "branchsite_classes_table_headers": [
        "Classe",        # largura 120 px
        "Proporção",     # largura 150 px
        "Background ω",  # largura 150 px
        "Foreground ω",  # largura 150 px
    ],

    # Valor exibido na coluna Background ω (requer parsing completo do arquivo)
    "branchsite_classes_bg_w_placeholder": "Ver arquivo",
    # Interpretação no rodapé da tabela de classes
    "branchsite_classes_footer": (
        "Valores de Foreground ω: valores elevados (> 1.0) indicam seleção positiva "
        "no ramo foreground para aquela classe de sítios"
    ),

    # ── Aba Branch Analysis (cladograma colorido por dN/dS) ──────────────
    "branch_tab_title":  "Branch Analysis — Cladograma por dN/dS",
    # Legenda de cores e instrução de interação
    "branch_tab_legend": (
        "Vermelho (w<1) · Amarelo (w=1) · Azul (w>1)"
        "   |   Clique num nó interno para girar"
    ),
    # Exibido quando não há dados do Branch Model no TSV
    "branch_no_data_title": "Sem dados do Branch Model",
    "branch_no_data_hint":  "Execute o modelo Branch com uma arvore etiquetada.",
    # Labels dos combos de gene e outgroup
    "branch_label_gene":     "Gene:",
    "branch_label_outgroup": "Outgroup:",
    # Valor padrão do combo de outgroup (nenhum re-enraizamento)
    "branch_outgroup_none":  "(nenhum)",
    # Botão de exportação do cladograma
    "branch_btn_export_png": "Exportar PNG",

    # ── Aba Export ────────────────────────────────────────────────────────
    # Título do cabeçalho da aba
    "export_tab_title": "Exportar Resultados",
    # Texto do botão de cada card de exportação
    "export_btn": "Exportar →",

    # Opções de exportação: lista de (título, descrição)
    # A ordem determina a ordem dos cards na interface.
    # NÃO altere a quantidade de itens — cada um corresponde a um método de exportação.
    "export_options": [
        ("Excel (.xlsx)",    "Tabela com formatação profissional"),   # → _export_excel
        ("CSV",              "Formato universal compatível"),          # → _export_csv
        ("Gráficos (PNG)",   "Exportar gráficos em alta resolução"),  # → _export_charts
        ("Relatório (HTML)", "Relatório completo interativo"),        # → _export_html
    ],
}
