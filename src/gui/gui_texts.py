"""Every interface text, in Portuguese (TEXTS_PT) and English (TEXTS_EN).
TEXTS reads from the active language; set_language() switches it."""

# ═══════════════════════════════════════════════════════════════════════════
# Portuguese
# ═══════════════════════════════════════════════════════════════════════════

TEXTS_PT: dict[str, object] = {

    # ── ModelConfigWindow ─────────────────────────────────────────
    "model_config_header":      "Modelo: {model_code}",
    "model_config_btn_cancel":  "Cancelar",
    "model_config_btn_save":    "Salvar",
    "model_status_default":     "Padrão",
    "model_status_configured":  "Configurado",

    # ── TreeLabelWindow ────────────────────────────────────────────
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

    # ── App (main window) ────────────────────────────────────
    "app_sidebar_title":    "EasyPAML",

    "section_files":   "ARQUIVOS",
    "section_results": "RESULTADOS",
    "section_config":  "CONFIGURAÇÕES",

    "btn_input_folder":  "Pasta de alinhamentos",
    "btn_tree_file":     "Arquivo de árvore",
    "btn_output_folder": "Pasta de Resultados",
    "label_not_selected": "Não selecionado",

    "btn_view_results":          "Ver Resultados",

    "label_omega_initial":    "ω (dN/dS) inicial:",
    "label_timeout":          "Tempo limite por modelo (min):",
    "label_timeout_auto":     "auto",
    "label_timeout_hint": (
        "Quanto tempo cada execução do codeml (um gene × um modelo) pode levar.\n\n"
        "Vazio = automático: o EasyPAML estima o tempo pelo modelo e pelo tamanho\n"
        "do gene (táxons × códons) e dá 5 vezes essa estimativa, no mínimo 30 min.\n"
        "Ex.: M8 com 30 táxons × 500 códons ≈ 33 min medidos → limite de 2,7 h.\n\n"
        "Preencha só se quiser um limite fixo (por exemplo, numa máquina muito\n"
        "lenta). Um codeml que pare de usar CPU por 5 min é encerrado de qualquer\n"
        "forma. Tabela das medições: docs/timing_benchmark.md."
    ),
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
    "label_omega_initial_hint": (
        "Valor de partida de ω (dN/dS) para a otimização por máxima verossimilhança\n"
        "(variável 'omega' do codeml; padrão do EasyPAML: 0,5). Não fixa ω: o codeml\n"
        "estima ω a partir dele (exceto nos nulos M8a e branch-site, em que ω = 1).\n\n"
        "O manual do PAML sugere repetir a análise com outros valores iniciais\n"
        "(ex.: 0,5 e 2) e conferir se o lnL chega ao mesmo valor: M7/M8 e\n"
        "branch-site podem ter problemas de convergência.\n\n"
        "Fonte: manual do PAML (pamlDOC: 'omega', 'Specifying initial values') e PAML FAQ."
    ),
    "label_ncatg_hint": (
        "Número de categorias usadas para discretizar a distribuição beta de ω\n"
        "nos modelos M7, M8 e M8a (variável 'ncatG' do codeml). Não afeta os\n"
        "outros modelos.\n\n"
        "Padrão: 10, o valor usado por Yang et al. (2000) para a beta e o que o\n"
        "codeml adota quando roda vários modelos num mesmo .ctl. Mais categorias\n"
        "aproximam melhor a distribuição e deixam a análise mais lenta.\n\n"
        "Fonte: manual do PAML (pamlDOC: 'NSsites' / 'ncatG')."
    ),
    "label_cpus_hint": (
        "Quantos genes o EasyPAML analisa ao mesmo tempo: um processo codeml\n"
        "por gene (o codeml usa um núcleo por processo). Os modelos de um mesmo\n"
        "gene rodam um depois do outro.\n\n"
        "Mais CPUs terminam o lote mais rápido, mas usam mais memória e deixam o\n"
        "computador mais lento para outras tarefas. Com um gene só, mais de 1 não\n"
        "acelera. Não muda os resultados."
    ),

    "label_ignore_stops":      "Ignorar stop codons",
    "label_ignore_stops_hint": (
        "Desligado (padrão): genes com stop codon no meio da\n"
        "sequência NÃO rodam e aparecem como FALHOU, com a\n"
        "sequência e a posição do códon.\n"
        "Ligado: o codeml roda e trata a coluna inteira do stop\n"
        "como dado ausente (com cleandata = 1 ela sai da análise).\n"
        "Um stop no último códon nunca impede a análise."
    ),

    "label_warm_start":      "Partir do M0 (warm start)",
    "label_warm_start_hint": (
        "Desligado (padrão): cada modelo é ajustado do zero.\n"
        "Ligado: o M0 é ajustado primeiro (mesmo sem estar marcado);\n"
        "seu κ e seus comprimentos de ramo viram valores iniciais\n"
        "dos outros modelos (fix_blength = 1), e cada um é ajustado\n"
        "a partir de ω inicial 0,2, 1,0 e 2,5, ficando o melhor lnL.\n"
        "Veja METHODS.md."
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
    "status_stops_template": "■ Stop codons: {n}",
    "msg_not_results_folder": "{path} não tem resultados do EasyPAML (analysis_summary.tsv). Escolha a pasta indicada em \u201cPasta de Resultados\u201d quando a análise rodou.",
    "label_neutral_models":  "Modelos nulos automáticos",

    "btn_run":           "Iniciar",
    "btn_pause":         "Pausar",
    "btn_resume":        "Retomar",
    "btn_stop":          "Parar",
    "btn_open_output":   "Abrir pasta de resultados",

    "log_header_title":    "LOG DE EXECUÇÃO",
    "log_header_subtitle": "uma mensagem por linha",

    "log_welcome": (
        "EasyPAML -- análise de seleção com PAML/codeml\n"
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

    # ── ResultsViewerWindow ───────────────────────────────────────
    "viewer_window_title":       "EasyPAML — Painel de Análise",
    "viewer_error_no_tsv":       "Arquivo analysis_summary.tsv não encontrado!",
    "viewer_error_run_analysis": "Execute uma análise para gerar resultados.",

    "viewer_header_title":    "EasyPAML  —  Resultados",
    "viewer_header_subtitle": "Análise de seleção  ·  codeml",
    "viewer_btn_recompute": "Recalcular resumos",

    "viewer_tab_lrt":                "LRT e p-valores",
    "viewer_tab_omega":              "ω > 1 Global",
    "viewer_tab_sites":              "Sítios sob seleção",
    "viewer_tab_branchsite_classes": "Classes branch-site",
    "viewer_tab_branch":             "Análise de ramos",

    "stats_total_genes":        "Total de Genes",
    "stats_models_run":         "Modelos Rodados",
    "stats_positive_selection": "Seleção Positiva Global",
    "stats_avg_omega":          "ω Médio",

    "lrt_no_comparisons": "Nenhuma comparação LRT disponível",
    "lrt_label_model":    "Teste:",


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


    "sites_label_model":    "Modelo:",
    "sites_label_gene":     "Gene:",
    "sites_label_analysis": "Análise:",
    "sites_label_filter":   "Filtrar Pr(w>1) ≥",
    "sites_btn_update":     "Atualizar",

    "sites_file_not_found": "Arquivo não encontrado: {filename}",
    "sites_parse_error":    "Erro ao ler o arquivo:\n{error}",
    "sites_no_sites":       "Nenhum sítio com Pr(ω>1) ≥ {threshold}",
    "sites_chart_title":    "● {n1} *  Pr(ω>1) ≥ 0,95      ◆ {n2} **  Pr(ω>1) ≥ 0,99",
    "sites_chart_removed":  "Hachurado: {n} códon(s) removido(s) pelo cleandata (gap, ambiguidade ou stop em alguma sequência); o codeml não os analisa.",
    "sites_show_all":       "Mostrar todos os genes",
    "sites_genes_significant": "{n} gene(s) significativo(s) (q < 0,05)",
    "sites_genes_shown":    "{n} gene(s)",
    "sites_no_significant_genes": "Nenhum gene com teste significativo (q < 0,05) para o {model}. Marque \"Mostrar todos os genes\" para ver os sítios dos outros.",

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

    "branch_no_data_title": "Sem dados do Branch Model",
    "branch_no_data_hint":  "Execute o modelo Branch com uma árvore marcada.",
    "branch_export": "Exportar árvore…",
    "branch_groups": "Grupos de ramos (modelo Branch)",
    "branch_col_group": "Grupo",
    "branch_col_n": "Ramos",
    "branch_col_w": "ω",
    "branch_background": "fundo",
    "branch_bs_info": "Branch-site, ramos marcados: {p2}% dos sítios na classe com ω > 1 (ω = {w2}); {n} sítio(s) com Pr(ω>1) ≥ 0,95 no BEB.",
    "branch_label_gene":     "Gene:",



    # ── TreeLabelWindow — mensagens de erro inline ───────────────────────
    "tree_err_no_biopython":  "[!] Biopython não instalado. Execute: pip install biopython",
    "tree_err_no_tree":       "Nenhuma árvore selecionada.",
    "tree_err_load":          "Erro ao carregar a árvore:\n{error}",

    # ── TreeLabelWindow dialogs ────────────────────────────────
    "tag_dialog_edit_title":  "Editar Tag",
    "tag_dialog_edit_prompt": "Ramo atual: {tag}\n\nDigite novo número ou 'remover':",
    "tag_dialog_new_title":   "Número da Tag",
    "tag_dialog_new_prompt":  "Digite o número da tag:\n(ex: 1 para #1, 2 para #2)",

    # ── App — CPU label ──────────────────────────────────────────────────
    "label_cpus_detected":    "(detectado: {n} núcleos)",

    # ── App — mensagens de log ───────────────────────────────────────────
    "log_no_tree_selected":   "Selecione um arquivo de árvore (.nwk) primeiro.\n",
    "log_no_output_folder":   "Selecione uma pasta de resultados primeiro.\n",
    
    "log_analysis_start":     "Iniciando a análise…\n",
    "log_analysis_stopped":   "Análise interrompida.\n",
    "log_analysis_done":      "Análise terminada.\n",

    # ── App file dialogs ────────────────────────────────────────
    "dialog_select_results_folder": "Escolha uma pasta de resultados do EasyPAML para abrir",

    # ── App — troca de idioma ────────────────────────────────────────────
    "lang_switch_message":  "O novo idioma é aplicado quando o EasyPAML reabre. As pastas e os modelos escolhidos são mantidos.",
    "lang_switch_confirm":  "Reabrir agora?",
    "msg_wait_for_run": "Espere a análise terminar (ou pare-a) antes de trocar o idioma ou o tema.",
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
    "msg_exported_to":          "Exportado para:\n{path}",
    "msg_export_err":           "Erro ao exportar:\n{error}",
    "msg_no_lrt":               "Nenhum resultado LRT encontrado.",
    "msg_excel_exported":       "Excel exportado com {n} aba(s):\n{path}",
    "msg_excel_err":            "Erro ao exportar Excel:\n{error}",
    "msg_html_exported":        "Relatório HTML exportado:\n{path}",
    "msg_html_err":             "Erro ao exportar HTML:\n{error}",

    # ── ResultsViewerWindow file dialogs ────────────────────────
    "dialog_save_as":           "Salvar como",

    # ── Charts (matplotlib) ─────────────────────────────────────────────
    "chart_whole_gene":         "gene inteiro",
    "chart_positive_class":     "classe positiva",
    "chart_hint":               "Passe o mouse sobre o eixo x para ver os genes. Laranja: q < 0,05.",

    # ── Novas chaves (0.3.0) ─────────────────────────────────────────────
    "picker_up": "Voltar para a pasta de cima (Backspace)",
    "picker_choose": "Escolher esta pasta",
    "picker_choose_named": "Escolher \u201c{name}\u201d",
    "picker_cancel": "Cancelar",
    "picker_new_folder": "Nova pasta",
    "picker_new_folder_prompt": "Nome da nova pasta:",
    "picker_file_name": "Nome:",
    "picker_open": "Abrir",
    "picker_save": "Salvar",
    "picker_overwrite": "\u201c{name}\u201d já existe. Substituir?",
    "picker_path_missing": "Não existe: {path}",
    "picker_hint_open": "Clique duas vezes numa pasta para abrir e num arquivo para escolhê-lo. Arquivos de outro tipo aparecem em cinza. Também dá para digitar o caminho acima e apertar Enter.",
    "picker_hint_save": "Escolha a pasta (clique duas vezes para abrir) e confira o nome do arquivo.",
    "picker_hint": "Clique duas vezes numa pasta para abrir; um clique a marca para escolher. Os arquivos aparecem em cinza só para você conferir o conteúdo. Também dá para digitar o caminho acima e apertar Enter.",
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
        "Métodos (parâmetros do codeml, LRT, correção BH): METHODS.md"
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
    "label_per_gene_trees": "{n} de {total} gene(s) com árvore própria (GENE.nwk) -- usada no lugar do arquivo de árvore",
    "label_output_created": "{name} (pasta criada)",
    "hint_select_files": "Escolha {missing} para ativar os modelos.",
    "hint_missing_parts": ("a pasta de alinhamentos", "o arquivo de árvore", "a pasta de resultados"),
    "hint_and": " e ",
    "progress_idle": "Nenhuma análise em andamento",
    "progress_template": "Gene {done} de {total}",
    "progress_done": "{done} de {total} genes processados",
    "progress_done_failed": "{ok} de {total} genes concluídos · {failed} falharam",
    "progress_running": "{pct}% · {done} de {total} genes prontos · {elapsed} decorridos{left} · rodando {what}",
    "progress_model_genes": "{model} em {n} genes",
    "progress_model_gene": "{model} em 1 gene",
    "progress_left": " · faltam cerca de {left}",
    "progress_stopped": "Parada por você · {ok} de {total} genes concluídos",
    "overwrite_title": "Pasta com resultados",
    "overwrite_text": (
        "Esta pasta já tem resultados. Genes e modelos com os mesmos dados e parâmetros são "
        "reaproveitados; os outros são rodados e substituem o que havia."
    ),
    "overwrite_yes": "Continuar",
    "overwrite_no": "Escolher outra pasta",
    "err_permission": "Você não tem permissão para gravar em {path}. Escolha uma pasta dentro da sua pasta pessoal.",
    "err_not_found": "{path} não existe (foi movido ou apagado?).",
    "err_disk_full": "Não há espaço livre no disco para gravar em {path}.",
    "msg_still_running": "A análise ainda está rodando. O painel de resultados abre sozinho quando ela terminar.",
    "stop_confirm_title": "Parar a análise?",
    "stop_confirm_text": (
        "Os modelos em execução serão interrompidos e os genes que ainda não terminaram "
        "ficam sem resultado. Os genes já concluídos são mantidos."
    ),
    "stop_confirm_yes": "Parar",
    "stop_confirm_no": "Continuar rodando",
    "chk_show_details": "Mostrar detalhes técnicos",
    "preflight_title": "Verificação dos dados",
    "preflight_ok": "Verificação dos dados: nenhum problema em {n} gene(s).",
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
    "preflight_btn_fix": "Voltar para corrigir os arquivos",
    "preflight_btn_continue": "Continuar mesmo assim",
    "preflight_continue_hint": (
        "Continuando, genes com erro aparecem como FALHOU e stop codons são tratados "
        "pelo codeml como dado ausente."
    ),
    "failures_title": "Genes que falharam",
    "failures_text": "{failed} de {total} gene(s) falharam:\n\n{items}\n\nDetalhes em {path}",
    "msg_output_folder_error": "Não foi possível criar a pasta de resultados:\n{error}",
    "dialog_choose_output": "Escolha a pasta para os resultados (ou crie uma com \"Nova pasta\")",
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
    "theme_label": "Tema",
    "theme_system": "Auto",
    "theme_light": "Claro",
    "theme_dark": "Escuro",
    "theme_restart": "O novo tema é aplicado quando o EasyPAML reabre. As pastas e os modelos escolhidos são mantidos. Reabrir agora?",
    "lang_restart_title": "Idioma",
    "viewer_tab_summary": "Resumo",
    "viewer_btn_open_output": "Abrir pasta de resultados",
    "stats_sig_genes": "Genes com seleção positiva (q < 0,05)",
    "loading_results": "Carregando os resultados…",
    "loading_app": "Carregando…",
    "stats_sig_short": "q < 0,05:",
    "stats_failed": "Genes que falharam",
    "summary_title": "Uma linha por gene e por teste",
    "summary_explain": (
        "q é o p corrigido para o número de genes (Benjamini-Hochberg); q < 0,05 é significativo. "
        "Com um gene só não há o que corrigir e q = p. "
        "\"ω médio\" é a média de todos os sítios no modelo alternativo do teste. Ele não decide o resultado: "
        "pode ficar abaixo de 1 com sítios sob seleção e acima de 1 sem um teste significativo. "
        "A evidência vem do q, da classe positiva e dos sítios."
    ),
    "summary_test": "Teste:",
    "summary_pick_gene": "Clique num gene para ver os detalhes; clique duas vezes para abrir os sítios.",
    "msg_open_folder_failed": "Não foi possível abrir um gerenciador de arquivos. Os resultados estão em:\n{path}\n(caminho copiado para a área de transferência)",
    "btn_try_example": "Testar com o exemplo",
    "log_example_loaded": "Exemplo carregado: 2 genes simulados com resposta conhecida (um com seleção positiva, outro sem; veja examples/quick/README.md) e o M8 ligado com seus modelos nulos. Clique em Run; leva poucos minutos. Resultados em {out}",
    "settings_defaults_note": "Valores padrão; dá para rodar sem mudar nada. Cada um é explicado no \"?\" ao lado.",
    "tree_choice_message": "Uma árvore para todos os genes, ou uma pasta com uma árvore por gene? Na pasta, cada árvore é pareada com o gene pelo nome do arquivo (GENE.nwk, GENE.fasta.treefile do IQ-TREE, RAxML_bestTree.GENE …). Dá para usar as duas: a árvore geral vale para os genes sem árvore própria.",
    "tree_choice_file": "Uma árvore para todos",
    "tree_choice_folder": "Pasta com uma árvore por gene",
    "label_tree_folder": "{folder}/ · {n} de {total} gene(s) pareado(s)",
    "label_tree_folder_general": "os outros usam {tree}",
    "pairing_title": "Árvores por gene",
    "pairing_summary": "{paired} de {total} gene(s) pareado(s) · {missing} sem árvore · {orphans} árvore(s) sem gene",
    "pairing_rules": "Os nomes são comparados sem extensões nem os prefixos do IQ-TREE e do RAxML (.treefile, .nwk, RAxML_bestTree.) e sem diferenciar maiúsculas. Um nome parecido nunca é pareado: renomeie o arquivo para ficar igual ao alinhamento e clique em Conferir de novo. O EasyPAML não renomeia seus arquivos.",
    "pairing_branch_note": "Branch e Branch-site: as marcas #1 precisam já estar em cada árvore por gene.",
    "pairing_col_gene": "Gene (alinhamento)",
    "pairing_col_tree": "Arquivo de árvore",
    "pairing_col_status": "Situação",
    "pairing_paired": "pareado",
    "pairing_no_tree": "nenhuma árvore com este nome",
    "pairing_did_you_mean": "você quis dizer {name}? renomeie o arquivo",
    "pairing_general_used": "· usa a árvore geral",
    "pairing_duplicate": "mais de uma árvore para este gene; deixe só uma",
    "pairing_orphan": "nenhum alinhamento com este nome",
    "pairing_open_folder": "Abrir a pasta de árvores",
    "pairing_check_again": "Conferir de novo",
    "pairing_continue": "Continuar",
    "answer_positive": "✓ Seleção positiva sustentada em {n} de {total} gene(s){genes}.",
    "answer_none": "Seleção positiva não sustentada em nenhum dos {total} gene(s).",
    "answer_neutral": "{n} significativo(s) só no M8 vs M7, o que sítios neutros (ω = 1) podem causar; o M8 vs M8a não confirma{genes}.",
    "answer_no_m8a": "{n} significativo(s) no M8 vs M7 sem o M8a para descartar sítios neutros{genes}.",
    "answer_weak": "{n} com teste significativo mas nenhum sítio com Pr(ω>1) ≥ 0,95 (sinal fraco){genes}.",
    "answer_sites": "{n} sítio(s)",
    "test_plain": {
        ("M8a", "M8"): "Significativo (q < 0,05): uma classe de sítios com ω > 1 explica os dados claramente melhor do que a mesma classe com ω = 1. É o teste mais rigoroso (Swanson et al. 2003).",
        ("M7", "M8"): "Significativo: o M8 explica os dados melhor que o M7, mas sítios neutros (ω = 1) sozinhos podem causar isso; confira o M8 vs M8a.",
        ("M1a", "M2a"): "Significativo: uma classe de sítios com ω > 1 explica os dados claramente melhor do que só ω < 1 e ω = 1.",
        ("Branch-site_null", "Branch-site"): "Significativo: há sítios com ω > 1 nos ramos marcados (foreground).",
        ("M0", "Branch"): "Significativo: o ω dos ramos marcados é diferente do resto da árvore; isso sozinho não mostra ω > 1.",
        ("M0", "M1a"): "Significativo: o ω varia entre os sítios. É um teste de variação, não de seleção positiva."},
    "chart_hide": "Ocultar gráfico",
    "chart_show": "Mostrar gráfico",
    "summary_warn_hint": "⚠ = o dado foi alterado antes da análise (clique no gene para ver).",
    "summary_open_sites": "Ver sítios",
    "summary_export_table": "Exportar tabela…",
    "summary_export_all": "Todos os testes (Excel)…",
    "summary_export_html": "Relatório HTML…",
    "summary_copy_genes": "Copiar genes significativos",
    "summary_genes_copied": "{n} gene(s) copiado(s), um por linha (para g:Profiler, PANTHER, DAVID…).",
    "summary_df_branch": "df = número de grupos de ramos marcados",
    "chart_kind_lrt": "Distribuição de 2Δℓ",
    "chart_kind_omega": "ω por gene",
    "chart_export": "Exportar gráfico…",
    "chart_mean_w": "ω médio ({model})",
    "col_gene": "Gene",
    "col_result": "Neste teste",
    "col_conclusion": "Evidência (todos os testes)",
    "col_mean_w": "ω médio ({model})",
    "col_w_pos": "ω classe +",
    "col_sites": "Sítios ≥0,95",
    "col_w_tags": "ω por marca",
    "result_failed": "falhou",
    "conclusion_short": {
        "positive": "Sustentada", "weak": "Sinal fraco", "neutral": "Não confirmada pelo M8a",
        "no_m8a": "Possível (sem M8a)", "none": "Não detectada", "failed": "Falhou", "sig": ""},
    "summary_col_hints": {
        "q": "p corrigido para o número de genes testados (Benjamini-Hochberg). q < 0,05 é significativo.",
        "p": "Probabilidade de uma melhora tão grande sem seleção positiva (distribuição χ²).",
        "lrt": "2Δℓ = 2 × (lnL alternativo − lnL nulo): quanto o modelo alternativo melhora o ajuste.",
        "mean_w": "Média de ω em todos os sítios no modelo alternativo. Não é critério de seleção positiva.",
        "w_pos": "ω da classe de sítios que pode estar sob seleção positiva.",
        "p1": "Proporção de sítios nessa classe.",
        "sites": "Sítios com Pr(ω>1) ≥ 0,95 no BEB.",
        "conclusion": "Evidência de seleção positiva juntando os testes de sítio do gene (M8 vs M8a, M8 vs M7, M2a vs M1a).",
        "result": "Significativo (q < 0,05) neste teste?",
        "lnl0": "Log-verossimilhança do modelo nulo; maior (menos negativo) é melhor ajuste.",
        "lnl1": "Log-verossimilhança do modelo alternativo.",
        "w_tags": "ω de cada grupo de ramos marcado (bg = ramos de fundo).",
        "df": "Graus de liberdade do teste: número de grupos marcados."},
    "test_hypotheses": {
        ("M8a", "M8"): "H₀ M8a: beta + classe com ω = 1 fixo   ·   H₁ M8: beta + classe com ω livre (Swanson et al. 2003)",
        ("M7", "M8"): "H₀ M7: beta entre 0 e 1   ·   H₁ M8: beta + classe com ω livre   ·   pode dar positivo só por sítios neutros",
        ("M1a", "M2a"): "H₀ M1a: purificação e neutralidade (ω ≤ 1)   ·   H₁ M2a: mais uma classe com ω > 1",
        ("Branch-site_null", "Branch-site"): "H₀: ω = 1 no foreground   ·   H₁: sítios com ω > 1 nos ramos marcados (#1)",
        ("M0", "Branch"): "H₀ M0: um ω para todos os ramos   ·   H₁ Branch: um ω por grupo de ramos marcado",
        ("M0", "M1a"): "H₀ M0: um ω para todos os sítios   ·   H₁ M1a: classes com ω < 1 e ω = 1"},
    "summary_sig": "{test}: significativo (p = {p}, q = {q}) — {sites} sítio(s) com Pr(ω>1) ≥ 0,95",
    "summary_nonsig": "{test}: não significativo (p = {p}, q = {q})",
    "summary_posclass": "classe positiva: ω = {w}, p₁ = {p1}",
    "conclusion_supported": "Seleção positiva sustentada: {tests} significativo(s) (q < 0,05){sites}.",
    "conclusion_sites": ", {n} sítio(s) com Pr(ω>1) ≥ 0,95",
    "conclusion_neutral": (
        "Não confirmada pelo M8a: M8 vs M7 é significativo mas M8a vs M8 não, então o sinal pode vir "
        "de sítios neutros (ω = 1)."
    ),
    "conclusion_none": "Seleção positiva não detectada (nenhum teste com q < 0,05).",
    "conclusion_weak": (
        "Sinal fraco: {tests} significativo(s), mas nenhum sítio com Pr(ω>1) ≥ 0,95{detail}. "
        "Sinais assim costumam vir de poucos códons mal alinhados; confira o alinhamento."
    ),
    "conclusion_weak_class": " (classe positiva com ω = {w} em {p1}% dos sítios)",
    "conclusion_no_m8a": (
        "Possível seleção positiva: só M8 vs M7 é significativo e o M8a não rodou, então sítios "
        "neutros (ω = 1) não foram descartados. Rode o M8a para conferir."
    ),
    "summary_failed": "FALHOU: {reason}",
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
    "sites_btn_figure": "Exportar figura…",
    "sites_copied": "{n} sítio(s) copiado(s) para a área de transferência.",
    "lrt_headers": ["Gene", "lnL nulo", "lnL alternativo", "2Δℓ", "p", "q (BH)", "ω classe + (p₁)", "Sig."],
    "lrt_plain": (
        "Cada linha compara dois modelos do mesmo gene. Um q pequeno (< 0,05) quer dizer que o "
        "modelo com seleção positiva explica os dados claramente melhor. Passe o mouse nos "
        "títulos das colunas para ver o que cada uma significa."
    ),
    "lrt_header_hints": [
        "Nome do gene (arquivo de alinhamento).",
        "lnL do modelo nulo (sem seleção positiva). lnL é o log da verossimilhança: quanto maior "
        "(menos negativo), melhor o modelo explica os dados.",
        "lnL do modelo alternativo (com uma classe que pode ter ω > 1).",
        "2Δℓ = 2 × (lnL alternativo − lnL nulo): quanto o modelo alternativo melhora o ajuste. "
        "É a estatística do teste (LRT).",
        "p: probabilidade de uma melhora assim aparecer sem seleção positiva, pela distribuição χ².",
        "q: o p corrigido para os muitos genes testados (Benjamini-Hochberg). É o número que "
        "decide: q < 0,05 é significativo.",
        "ω e proporção (p₁) da classe de sítios que pode estar sob seleção positiva.",
        "Significativo (q < 0,05)?",
    ],
    "lrt_sig_yes": "sim",
    "lrt_sig_no": "não",
    "summary_verdict_sig": "significativo",
    "summary_verdict_nonsig": "não significativo",
    "summary_verdict_failed": "FALHOU",
    "summary_verdict_warning": "AVISO",
    "summary_sites_n": "{n} sítio(s) com Pr(ω>1) ≥ 0,95",
    # ── Main window: steps, cards, summary ──
    "app_main_title":     "Análise de seleção",
    "app_main_subtitle":  "Processamento do codeml em lote",
    "step_data":          "Dados",
    "step_models":        "Modelos",
    "step_settings":      "Configurações avançadas",
    "step_models_none":   "Ligue os modelos ao lado",
    "step_models_count":  "{n} selecionado(s)",
    "slot_tree_per_gene": "{n} árvore(s) por gene (GENE.nwk), não precisa escolher outra",
    "slot_choose":        "Escolher",
    "slot_change":        "Trocar",
    "run_summary_no_data":   "Escolha os dados na etapa 1",
    "run_summary_no_models": "nenhum modelo ligado",
    "run_summary_genes":     "{n} gene(s)",
    "run_summary_cpus":      "{n} CPU(s)",
    "settings_show":       "Mostrar",
    "settings_hide":       "Ocultar",
    "log_collapse":        "Recolher",
    "log_expand":          "Expandir",
    "run_summary_models":  "{n} modelo(s)",
    "run_tests_label":     "Testes LRT: {tests}",
    "run_tests_none":      ("Nenhum teste LRT será feito: cada teste precisa de um par nulo + "
                            "alternativo (ex.: M1a e M2a, M7 e M8)."),
    "tile_auto":           "Entra automaticamente como nulo do {alt}. Clique para deixá-lo de fora.",
    "tile_auto_off":       "Deixado de fora: o teste {alt} vs {null} não roda. Clique para incluí-lo.",
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

    "section_files":   "FILES",
    "section_results": "RESULTS",
    "section_config":  "SETTINGS",

    "btn_input_folder":   "Alignments folder",
    "btn_tree_file":      "Tree file",
    "btn_output_folder":  "Output folder",
    "label_not_selected": "Not selected",

    "btn_view_results":          "View Results",

    "label_omega_initial":    "Initial ω (dN/dS):",
    "label_timeout":          "Time limit per model (min):",
    "label_timeout_auto":     "auto",
    "label_timeout_hint": (
        "How long each codeml run (one gene × one model) may take.\n\n"
        "Empty = automatic: EasyPAML estimates the time from the model and the\n"
        "gene size (taxa × codons) and allows 5 times that, at least 30 min.\n"
        "E.g. M8 with 30 taxa × 500 codons ≈ 33 min measured → 2.7 h limit.\n\n"
        "Fill it in only if you want a fixed limit (e.g. on a very slow\n"
        "machine). A codeml that stops using CPU for 5 min is stopped anyway.\n"
        "Measurements: docs/timing_benchmark.md."
    ),
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
    "label_omega_initial_hint": (
        "Starting value of ω (dN/dS) for the maximum-likelihood optimisation\n"
        "(codeml variable 'omega'; EasyPAML default: 0.5). It does not fix ω:\n"
        "codeml estimates ω from it (except in the M8a and branch-site nulls,\n"
        "where ω = 1).\n\n"
        "The PAML manual suggests re-running with other initial values (e.g. 0.5\n"
        "and 2) and checking that lnL reaches the same value: M7/M8 and\n"
        "branch-site models can have convergence problems.\n\n"
        "Source: PAML manual (pamlDOC: 'omega', 'Specifying initial values') and PAML FAQ."
    ),
    "label_ncatg_hint": (
        "Number of categories used to discretise the beta distribution of ω in\n"
        "models M7, M8 and M8a (codeml variable 'ncatG'). It does not affect the\n"
        "other models.\n\n"
        "Default: 10, the value used by Yang et al. (2000) for the beta and the one\n"
        "codeml uses when several models run from one .ctl. More categories\n"
        "approximate the distribution better and make the analysis slower.\n\n"
        "Source: PAML manual (pamlDOC: 'NSsites' / 'ncatG')."
    ),
    "label_cpus_hint": (
        "How many genes EasyPAML analyses at the same time: one codeml process\n"
        "per gene (codeml uses one core per process). The models of a gene run\n"
        "one after the other.\n\n"
        "More CPUs finish the batch sooner but use more memory and slow the\n"
        "computer down for other tasks. With a single gene, more than 1 does not\n"
        "help. It does not change the results."
    ),

    "label_ignore_stops":      "Ignore stop codons",
    "label_ignore_stops_hint": (
        "Off (default): genes with a stop codon inside the\n"
        "sequence do NOT run and are shown as FAILED, with the\n"
        "sequence and codon position.\n"
        "On: codeml runs and treats the whole stop column as\n"
        "missing data (with cleandata = 1 it is removed).\n"
        "A stop in the last codon never blocks the analysis."
    ),

    "label_warm_start":      "Start from M0 (warm start)",
    "label_warm_start_hint": (
        "Off (default): every model is fitted from scratch.\n"
        "On: M0 is fitted first (even if it is not selected);\n"
        "its κ and branch lengths become the starting values\n"
        "of the other models (fix_blength = 1), and each one is\n"
        "fitted from initial ω 0.2, 1.0 and 2.5, keeping the best\n"
        "lnL. See METHODS.md."
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
    "status_stops_template": "■ Stop codons: {n}",
    "msg_not_results_folder": "{path} has no EasyPAML results (analysis_summary.tsv). Choose the folder given as \u201cOutput folder\u201d when the analysis was run.",
    "label_neutral_models":  "Automatic null models",

    "btn_run":           "Run",
    "btn_pause":         "Pause",
    "btn_resume":        "Resume",
    "btn_stop":          "Stop",
    "btn_open_output":   "Open results folder",

    "log_header_title":    "EXECUTION LOG",
    "log_header_subtitle": "one message per line",

    "log_welcome": (
        "EasyPAML -- selection analysis with PAML/codeml\n"
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
    "viewer_header_subtitle": "Selection Analysis  ·  codeml",
    "viewer_btn_recompute": "Recompute summaries",

    "viewer_tab_lrt":                "LRT and p-values",
    "viewer_tab_omega":              "ω > 1 Global",
    "viewer_tab_sites":              "Positive Sites",
    "viewer_tab_branchsite_classes": "Branch-site Classes",
    "viewer_tab_branch":             "Branch Analysis",

    "stats_total_genes":        "Total Genes",
    "stats_models_run":         "Models Run",
    "stats_positive_selection": "Global Positive Selection",
    "stats_avg_omega":          "Average ω",

    "lrt_no_comparisons": "No LRT comparison available",
    "lrt_label_model":    "Test:",


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


    "sites_label_model":    "Model:",
    "sites_label_gene":     "Gene:",
    "sites_label_analysis": "Analysis:",
    "sites_label_filter":   "Filter Pr(w>1) ≥",
    "sites_btn_update":     "Update",

    "sites_file_not_found": "File not found: {filename}",
    "sites_parse_error":    "Error reading the file:\n{error}",
    "sites_no_sites":       "No site with Pr(ω>1) ≥ {threshold}",
    "sites_chart_title":    "● {n1} *  Pr(ω>1) ≥ 0.95      ◆ {n2} **  Pr(ω>1) ≥ 0.99",
    "sites_chart_removed":  "Hatched: {n} codon(s) removed by cleandata (gap, ambiguity or stop codon in some sequence); codeml does not analyse them.",
    "sites_show_all":       "Show all genes",
    "sites_genes_significant": "{n} significant gene(s) (q < 0.05)",
    "sites_genes_shown":    "{n} gene(s)",
    "sites_no_significant_genes": "No gene with a significant test (q < 0.05) for {model}. Tick \"Show all genes\" to see the sites of the others.",

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

    "branch_no_data_title": "No Branch Model data",
    "branch_no_data_hint":  "Run the Branch model with a labeled tree.",
    "branch_export": "Export tree…",
    "branch_groups": "Branch groups (Branch model)",
    "branch_col_group": "Group",
    "branch_col_n": "Branches",
    "branch_col_w": "ω",
    "branch_background": "background",
    "branch_bs_info": "Branch-site, labelled branches: {p2}% of sites in the class with ω > 1 (ω = {w2}); {n} site(s) with Pr(ω>1) ≥ 0.95 in BEB.",
    "branch_label_gene":     "Gene:",



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
    
    "log_analysis_start":     "Starting the analysis…\n",
    "log_analysis_stopped":   "Analysis stopped.\n",
    "log_analysis_done":      "Analysis finished.\n",

    # ── App — file dialogs ───────────────────────────────────────────────
    "dialog_select_results_folder": "Choose an EasyPAML results folder to open",

    # ── App — language switch ────────────────────────────────────────────
    "lang_switch_message":  "The new language is applied when EasyPAML reopens. The folders and models you chose are kept.",
    "lang_switch_confirm":  "Reopen now?",
    "msg_wait_for_run": "Wait for the analysis to finish (or stop it) before changing the language or theme.",
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
    "msg_exported_to":          "Exported to:\n{path}",
    "msg_export_err":           "Export error:\n{error}",
    "msg_no_lrt":               "No LRT results found.",
    "msg_excel_exported":       "Excel exported with {n} sheet(s):\n{path}",
    "msg_excel_err":            "Error exporting Excel:\n{error}",
    "msg_html_exported":        "HTML report exported:\n{path}",
    "msg_html_err":             "Error exporting HTML:\n{error}",

    # ── ResultsViewerWindow — file dialogs ───────────────────────────────
    "dialog_save_as":           "Save as",

    # ── Charts (matplotlib) ───────────────────────────────────────────────
    "chart_whole_gene":         "whole gene",
    "chart_positive_class":     "positive class",
    "chart_hint":               "Move the mouse along the x axis to see the genes. Orange: q < 0.05.",

    # ── New keys (0.3.0) ─────────────────────────────────────────────────
    "picker_up": "Back to the parent folder (Backspace)",
    "picker_choose": "Choose this folder",
    "picker_choose_named": "Choose \u201c{name}\u201d",
    "picker_cancel": "Cancel",
    "picker_new_folder": "New folder",
    "picker_new_folder_prompt": "Name of the new folder:",
    "picker_file_name": "Name:",
    "picker_open": "Open",
    "picker_save": "Save",
    "picker_overwrite": "\u201c{name}\u201d already exists. Replace it?",
    "picker_path_missing": "Not found: {path}",
    "picker_hint_open": "Double-click a folder to open it and a file to choose it. Files of other types are shown in grey. You can also type a path above and press Enter.",
    "picker_hint_save": "Choose the folder (double-click to open it) and check the file name.",
    "picker_hint": "Double-click a folder to open it; one click marks it to choose. Files are shown in grey only so you can check the contents. You can also type a path above and press Enter.",
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
        "Methods (codeml parameters, LRT, BH correction): METHODS.md"
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
    "label_per_gene_trees": "{n} of {total} gene(s) with their own tree (GENE.nwk) -- used instead of the tree file",
    "label_output_created": "{name} (folder created)",
    "hint_select_files": "Choose {missing} to enable the models.",
    "hint_missing_parts": ("the alignments folder", "the tree file", "the results folder"),
    "hint_and": " and ",
    "progress_idle": "No analysis running",
    "progress_template": "Gene {done} of {total}",
    "progress_done": "{done} of {total} genes processed",
    "progress_done_failed": "{ok} of {total} genes completed · {failed} failed",
    "progress_running": "{pct}% · {done} of {total} genes done · {elapsed} elapsed{left} · running {what}",
    "progress_model_genes": "{model} on {n} genes",
    "progress_model_gene": "{model} on 1 gene",
    "progress_left": " · about {left} left",
    "progress_stopped": "Stopped by you · {ok} of {total} genes completed",
    "overwrite_title": "Folder with results",
    "overwrite_text": (
        "This folder already has results. Genes and models with the same data and settings "
        "are reused; the others run and replace what was there."
    ),
    "overwrite_yes": "Continue",
    "overwrite_no": "Choose another folder",
    "err_permission": "You do not have permission to write in {path}. Choose a folder inside your home folder.",
    "err_not_found": "{path} does not exist (was it moved or deleted?).",
    "err_disk_full": "There is no free disk space to write in {path}.",
    "msg_still_running": "The analysis is still running. The results panel opens by itself when it ends.",
    "stop_confirm_title": "Stop the analysis?",
    "stop_confirm_text": (
        "The models that are running will be interrupted and genes that have not finished "
        "will have no results. Genes already completed are kept."
    ),
    "stop_confirm_yes": "Stop",
    "stop_confirm_no": "Keep running",
    "chk_show_details": "Show technical details",
    "preflight_title": "Data check",
    "preflight_ok": "Data check: no problems in {n} gene(s).",
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
    "preflight_btn_fix": "Go back to fix the files",
    "preflight_btn_continue": "Continue anyway",
    "preflight_continue_hint": (
        "If you continue, genes with errors are shown as FAILED and stop codons are "
        "treated by codeml as missing data."
    ),
    "failures_title": "Genes that failed",
    "failures_text": "{failed} of {total} gene(s) failed:\n\n{items}\n\nDetails in {path}",
    "msg_output_folder_error": "Could not create the results folder:\n{error}",
    "dialog_choose_output": "Choose the folder for the results (or create one with \"New folder\")",
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
    "theme_label": "Theme",
    "theme_system": "Auto",
    "theme_light": "Light",
    "theme_dark": "Dark",
    "theme_restart": "The new theme is applied when EasyPAML reopens. The folders and models you chose are kept. Reopen now?",
    "lang_restart_title": "Language",
    "viewer_tab_summary": "Summary",
    "viewer_btn_open_output": "Open results folder",
    "stats_sig_genes": "Genes with positive selection (q < 0.05)",
    "loading_results": "Loading results…",
    "loading_app": "Loading…",
    "stats_sig_short": "q < 0.05:",
    "stats_failed": "Genes that failed",
    "summary_title": "One line per gene and test",
    "summary_explain": (
        "q is p corrected for the number of genes (Benjamini-Hochberg); q < 0.05 is significant. "
        "With a single gene there is nothing to correct and q = p. "
        "\"Mean ω\" is the average over all sites under the test's alternative model. It does not decide the result: "
        "it can stay below 1 with sites under selection and exceed 1 without a significant test. "
        "The evidence comes from q, the positive class and the sites."
    ),
    "summary_test": "Test:",
    "summary_pick_gene": "Click a gene to see the details; double-click to open its sites.",
    "msg_open_folder_failed": "Could not open a file manager. The results are in:\n{path}\n(path copied to the clipboard)",
    "btn_try_example": "Try the example",
    "log_example_loaded": "Example loaded: 2 simulated genes with a known answer (one under positive selection, one not; see examples/quick/README.md) and M8 on with its null models. Click Run; it takes a few minutes. Results in {out}",
    "settings_defaults_note": "Default values; you can run without changing anything. Each one is explained by its \"?\".",
    "tree_choice_message": "One tree for all genes, or a folder with one tree per gene? In the folder, each tree is paired with its gene by file name (GENE.nwk, GENE.fasta.treefile from IQ-TREE, RAxML_bestTree.GENE …). Both can be used: the general tree is for the genes without a tree of their own.",
    "tree_choice_file": "One tree for all genes",
    "tree_choice_folder": "Folder with one tree per gene",
    "label_tree_folder": "{folder}/ · {n} of {total} gene(s) paired",
    "label_tree_folder_general": "the others use {tree}",
    "pairing_title": "Trees per gene",
    "pairing_summary": "{paired} of {total} gene(s) paired · {missing} without a tree · {orphans} tree(s) without a gene",
    "pairing_rules": "Names are compared without extensions and the IQ-TREE and RAxML prefixes (.treefile, .nwk, RAxML_bestTree.) and ignoring case. A near name is never paired: rename the file to match the alignment and click Check again. EasyPAML does not rename your files.",
    "pairing_branch_note": "Branch and Branch-site: the #1 labels must already be in each per-gene tree.",
    "pairing_col_gene": "Gene (alignment)",
    "pairing_col_tree": "Tree file",
    "pairing_col_status": "Status",
    "pairing_paired": "paired",
    "pairing_no_tree": "no tree with this name",
    "pairing_did_you_mean": "did you mean {name}? rename the file",
    "pairing_general_used": "· general tree used",
    "pairing_duplicate": "more than one tree for this gene; keep only one",
    "pairing_orphan": "no alignment with this name",
    "pairing_open_folder": "Open trees folder",
    "pairing_check_again": "Check again",
    "pairing_continue": "Continue",
    "answer_positive": "✓ Positive selection supported in {n} of {total} gene(s){genes}.",
    "answer_none": "Positive selection supported in none of the {total} gene(s).",
    "answer_neutral": "{n} significant only in M8 vs M7, which neutral sites (ω = 1) can cause; M8 vs M8a does not confirm{genes}.",
    "answer_no_m8a": "{n} significant in M8 vs M7 without M8a to rule out neutral sites{genes}.",
    "answer_weak": "{n} with a significant test but no site at Pr(ω>1) ≥ 0.95 (weak signal){genes}.",
    "answer_sites": "{n} site(s)",
    "test_plain": {
        ("M8a", "M8"): "Significant (q < 0.05): a class of sites with ω > 1 explains the data clearly better than the same class with ω = 1. This is the stricter test (Swanson et al. 2003).",
        ("M7", "M8"): "Significant: M8 explains the data better than M7, but neutral sites (ω = 1) alone can cause this; check M8 vs M8a.",
        ("M1a", "M2a"): "Significant: a class of sites with ω > 1 explains the data clearly better than ω < 1 and ω = 1 only.",
        ("Branch-site_null", "Branch-site"): "Significant: some sites have ω > 1 on the labelled (foreground) branches.",
        ("M0", "Branch"): "Significant: ω on the labelled branches differs from the rest of the tree; this alone does not show ω > 1.",
        ("M0", "M1a"): "Significant: ω varies among sites. This tests variation, not positive selection."},
    "chart_hide": "Hide chart",
    "chart_show": "Show chart",
    "summary_warn_hint": "⚠ = the data was changed before the run (click the gene to see how).",
    "summary_open_sites": "View sites",
    "summary_export_table": "Export table…",
    "summary_export_all": "All tests (Excel)…",
    "summary_export_html": "HTML report…",
    "summary_copy_genes": "Copy significant genes",
    "summary_genes_copied": "{n} gene(s) copied, one per line (for g:Profiler, PANTHER, DAVID…).",
    "summary_df_branch": "df = number of labelled branch groups",
    "chart_kind_lrt": "2Δℓ distribution",
    "chart_kind_omega": "ω per gene",
    "chart_export": "Export chart…",
    "chart_mean_w": "mean ω ({model})",
    "col_gene": "Gene",
    "col_result": "This test",
    "col_conclusion": "Evidence (all tests)",
    "col_mean_w": "mean ω ({model})",
    "col_w_pos": "ω pos. class",
    "col_sites": "Sites ≥0.95",
    "col_w_tags": "ω per label",
    "result_failed": "failed",
    "conclusion_short": {
        "positive": "Supported", "weak": "Weak signal", "neutral": "Not confirmed by M8a",
        "no_m8a": "Possible (no M8a)", "none": "Not detected", "failed": "Failed", "sig": ""},
    "summary_col_hints": {
        "q": "p corrected for the number of genes tested (Benjamini-Hochberg). q < 0.05 is significant.",
        "p": "The chance of an improvement this large without positive selection (χ² distribution).",
        "lrt": "2Δℓ = 2 × (alternative lnL − null lnL): how much the alternative model improves the fit.",
        "mean_w": "Average ω over all sites under the alternative model. Not a positive selection criterion.",
        "w_pos": "ω of the class of sites that may be under positive selection.",
        "p1": "Proportion of sites in that class.",
        "sites": "Sites with Pr(ω>1) ≥ 0.95 in BEB.",
        "conclusion": "Evidence for positive selection from all site tests of the gene (M8 vs M8a, M8 vs M7, M2a vs M1a).",
        "result": "Significant (q < 0.05) in this test?",
        "lnl0": "Log-likelihood of the null model; higher (less negative) is a better fit.",
        "lnl1": "Log-likelihood of the alternative model.",
        "w_tags": "ω of each labelled branch group (bg = background branches).",
        "df": "Degrees of freedom of the test: the number of labelled groups."},
    "test_hypotheses": {
        ("M8a", "M8"): "H₀ M8a: beta + a class with ω = 1 fixed   ·   H₁ M8: beta + a class with free ω (Swanson et al. 2003)",
        ("M7", "M8"): "H₀ M7: beta between 0 and 1   ·   H₁ M8: beta + a class with free ω   ·   can be positive from neutral sites alone",
        ("M1a", "M2a"): "H₀ M1a: purifying and neutral (ω ≤ 1)   ·   H₁ M2a: an extra class with ω > 1",
        ("Branch-site_null", "Branch-site"): "H₀: ω = 1 on the foreground   ·   H₁: sites with ω > 1 on the labelled branches (#1)",
        ("M0", "Branch"): "H₀ M0: one ω for all branches   ·   H₁ Branch: one ω per labelled branch group",
        ("M0", "M1a"): "H₀ M0: one ω for all sites   ·   H₁ M1a: classes with ω < 1 and ω = 1"},
    "summary_sig": "{test}: significant (p = {p}, q = {q}) — {sites} site(s) with Pr(ω>1) ≥ 0.95",
    "summary_nonsig": "{test}: not significant (p = {p}, q = {q})",
    "summary_posclass": "positive class: ω = {w}, p₁ = {p1}",
    "conclusion_supported": "Positive selection supported: {tests} significant (q < 0.05){sites}.",
    "conclusion_sites": ", {n} site(s) with Pr(ω>1) ≥ 0.95",
    "conclusion_neutral": (
        "Not confirmed by M8a: M8 vs M7 is significant but M8a vs M8 is not, so the signal may "
        "come from neutral sites (ω = 1)."
    ),
    "conclusion_none": "Positive selection not detected (no test with q < 0.05).",
    "conclusion_weak": (
        "Weak signal: {tests} significant, but no site has Pr(ω>1) ≥ 0.95{detail}. Signals "
        "like this often come from a few misaligned codons; check the alignment."
    ),
    "conclusion_weak_class": " (positive class with ω = {w} on {p1}% of sites)",
    "conclusion_no_m8a": (
        "Possible positive selection: only M8 vs M7 is significant and M8a was not run, so "
        "neutral sites (ω = 1) are not ruled out. Run M8a to check."
    ),
    "summary_failed": "FAILED: {reason}",
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
    "sites_btn_figure": "Export figure…",
    "sites_copied": "{n} site(s) copied to the clipboard.",
    "lrt_headers": ["Gene", "lnL null", "lnL alternative", "2Δℓ", "p", "q (BH)", "positive-class ω (p₁)", "Sig."],
    "lrt_plain": (
        "Each line compares two models of the same gene. A small q (< 0.05) means the model "
        "with positive selection explains the data clearly better. Move the mouse over a column "
        "title to see what it means."
    ),
    "lrt_header_hints": [
        "Gene name (alignment file).",
        "lnL of the null model (no positive selection). lnL is the log-likelihood: the higher "
        "(less negative), the better the model explains the data.",
        "lnL of the alternative model (with a class of sites that may have ω > 1).",
        "2Δℓ = 2 × (alternative lnL − null lnL): how much the alternative model improves the "
        "fit. It is the test statistic (LRT).",
        "p: the chance of an improvement this large without positive selection, from the χ² "
        "distribution.",
        "q: p corrected for the many genes tested (Benjamini-Hochberg). This is the number that "
        "decides: q < 0.05 is significant.",
        "ω and proportion (p₁) of the class of sites that may be under positive selection.",
        "Significant (q < 0.05)?",
    ],
    "lrt_sig_yes": "yes",
    "lrt_sig_no": "no",
    "summary_verdict_sig": "significant",
    "summary_verdict_nonsig": "not significant",
    "summary_verdict_failed": "FAILED",
    "summary_verdict_warning": "WARNING",
    "summary_sites_n": "{n} site(s) with Pr(ω>1) ≥ 0.95",
    # ── Main window (2nd visual pass): steps, tiles, summary ──
    "app_main_title":     "Selection Analysis",
    "app_main_subtitle":  "codeml batch processing",
    "step_data":          "Data",
    "step_models":        "Models",
    "step_settings":      "Advanced settings",
    "step_models_none":   "Switch models on at the right",
    "step_models_count":  "{n} selected",
    "slot_tree_per_gene": "{n} per-gene tree(s) (GENE.nwk), no need to choose one",
    "slot_choose":        "Choose",
    "slot_change":        "Change",
    "run_summary_no_data":   "Choose the data in step 1",
    "run_summary_no_models": "no model switched on",
    "run_summary_genes":     "{n} gene(s)",
    "run_summary_cpus":      "{n} CPU(s)",
    "settings_show":       "Show",
    "settings_hide":       "Hide",
    "log_collapse":        "Collapse",
    "log_expand":          "Expand",
    "run_summary_models":  "{n} model(s)",
    "run_tests_label":     "LRT tests: {tests}",
    "run_tests_none":      ("No LRT test will be run: each test needs a null + alternative pair "
                            "(e.g. M1a and M2a, M7 and M8)."),
    "tile_auto":           "Included automatically as the null of {alt}. Click to leave it out.",
    "tile_auto_off":       "Left out: the {alt} vs {null} test will not run. Click to include it.",
}


# ═══════════════════════════════════════════════════════════════════════════
# TEXTS["key"] reads from the active language
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
    """Short text in the active language, for strings built in code."""
    return pt if _current_lang == 'pt' else en


def get_language() -> str:
    """Active language code ('pt' or 'en')."""
    return _current_lang


class _TextProxy:
    """TEXTS['key'] always reads from the active language."""

    def __getitem__(self, key: str):
        return _AVAILABLE.get(_current_lang, TEXTS_PT)[key]

    def get(self, key: str, default=None):
        return _AVAILABLE.get(_current_lang, TEXTS_PT).get(key, default)

    def __contains__(self, key: str) -> bool:
        return key in _AVAILABLE.get(_current_lang, TEXTS_PT)


TEXTS = _TextProxy()
