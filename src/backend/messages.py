"""
Mensagens que o backend mostra ao usuário (log da janela e do CLI), em
português e inglês. A GUI chama set_language() com o idioma da interface; o
CLI usa --lang (padrão: inglês, ou português se o sistema estiver em pt).
"""

import locale
import os

_LANG = 'en'

MESSAGES = {
    'run_start': {
        'pt': "{n} gene(s) × {m} modelo(s) · {w} processo(s) em paralelo · codeml {version}",
        'en': "{n} gene(s) × {m} model(s) · {w} parallel worker(s) · codeml {version}",
    },
    'run_models': {
        'pt': "Modelos: {models} · CodonFreq {codonfreq} · ncatG {ncatg} · cleandata {cleandata}",
        'en': "Models: {models} · CodonFreq {codonfreq} · ncatG {ncatg} · cleandata {cleandata}",
    },
    'no_codeml': {
        'pt': "codeml não encontrado. Linux: sudo apt install paml · Windows: bin\\codeml.exe · "
              "ou defina EASYPAML_CODEML com o caminho do executável.",
        'en': "codeml not found. Linux: sudo apt install paml · Windows: bin\\codeml.exe · "
              "or set EASYPAML_CODEML to the executable path.",
    },
    'duplicate_ignored': {
        'pt': "{ignored} ignorado: é o mesmo gene que {used}.",
        'en': "{ignored} ignored: same gene as {used}.",
    },
    'gene_start': {
        'pt': "[{i}/{n}] {gene}",
        'en': "[{i}/{n}] {gene}",
    },
    'model_running': {
        'pt': "    {gene} · {model}: rodando…",
        'en': "    {gene} · {model}: running…",
    },
    'model_ok': {
        'pt': "    {gene} · {model}: lnL = {lnl:.4f}  ({t:.1f} s)",
        'en': "    {gene} · {model}: lnL = {lnl:.4f}  ({t:.1f} s)",
    },
    'model_failed': {
        'pt': "FALHOU: {gene} · {model}: {reason}",
        'en': "FAILED: {gene} · {model}: {reason}",
    },
    'gene_failed': {
        'pt': "FALHOU: {gene}: {reason}",
        'en': "FAILED: {gene}: {reason}",
    },
    'gene_ok': {
        'pt': "[OK] {gene}",
        'en': "[OK] {gene}",
    },
    'warn_stops_masked': {
        'pt': "AVISO: {gene}: {count} stop codon(s) ({details}); o codeml trata essas colunas como dado ausente.",
        'en': "WARNING: {gene}: {count} stop codon(s) ({details}); codeml treats those columns as missing data.",
    },
    'warn_excluded': {
        'pt': "AVISO: {gene}: {names} não está(ão) na árvore e foi(ram) EXCLUÍDA(S) da análise.",
        'en': "WARNING: {gene}: {names} not in the tree and EXCLUDED from the analysis.",
    },
    'warn_pruned': {
        'pt': "{gene}: {count} táxon(s) da árvore ausentes neste gene foram podados.",
        'en': "{gene}: {count} tree taxon/taxa absent from this gene were pruned.",
    },
    'warn_sitemap': {
        'pt': "AVISO: {gene} · {model}: o codeml analisou {codeml} códons e o EasyPAML esperava {expected}; "
              "a numeração original dos sítios não será mostrada.",
        'en': "WARNING: {gene} · {model}: codeml analysed {codeml} codons but EasyPAML expected {expected}; "
              "original site numbering will not be shown.",
    },
    'reason_timeout': {
        'pt': "o codeml passou do tempo limite ({s} s); para aumentar, use Configurações > "
              "Tempo limite ou --timeout",
        'en': "codeml exceeded the time limit ({s} s); to raise it, use Settings > Time limit "
              "or --timeout",
    },
    'reason_idle': {
        'pt': "o codeml ficou {s} s sem usar CPU (provavelmente esperando uma resposta); última linha: {line}",
        'en': "codeml used no CPU for {s} s (probably waiting for input); last line: {line}",
    },
    'reason_rc': {
        'pt': "o codeml terminou com código {rc}: {line}",
        'en': "codeml exited with code {rc}: {line}",
    },
    'reason_no_lnl': {
        'pt': "o codeml não informou a verossimilhança (lnL): {line}",
        'en': "codeml did not report a likelihood (lnL): {line}",
    },
    'reason_no_output': {
        'pt': "o codeml não gerou o arquivo de saída: {line}",
        'en': "codeml produced no output file: {line}",
    },
    'reason_stop_codons': {
        'pt': "stop codon(s) no meio da sequência: {details}. Corrija o alinhamento ou ative "
              "'Ignorar stop codons' para o codeml tratar essas colunas como dado ausente.",
        'en': "internal stop codon(s): {details}. Fix the alignment or enable 'Ignore stop codons' "
              "so codeml treats those columns as missing data.",
    },
    'reason_stopped': {
        'pt': "interrompido pelo usuário",
        'en': "stopped by the user",
    },
    'reason_exception': {
        'pt': "erro interno: {error}",
        'en': "internal error: {error}",
    },
    'reason_no_codeml': {
        'pt': "codeml não encontrado",
        'en': "codeml not found",
    },
    'reason_models': {
        'pt': "modelo(s) que falharam: {models}",
        'en': "failed model(s): {models}",
    },
    'lrt_start': {
        'pt': "Calculando testes de razão de verossimilhança (LRT)…",
        'en': "Computing likelihood ratio tests (LRT)…",
    },
    'lrt_none': {
        'pt': "Nenhum par de modelos aninhados para o LRT.",
        'en': "No pair of nested models for the LRT.",
    },
    'lrt_pair_done': {
        'pt': "LRT {null} vs {alt}: {n} gene(s), {sig} significativo(s) (q < 0,05)",
        'en': "LRT {null} vs {alt}: {n} gene(s), {sig} significant (q < 0.05)",
    },
    'summary_ok': {
        'pt': "ANÁLISE CONCLUÍDA: {ok} de {n} genes concluídos ({minutes:.1f} min)",
        'en': "ANALYSIS COMPLETE: {ok} of {n} genes completed ({minutes:.1f} min)",
    },
    'summary_failed': {
        'pt': "ANÁLISE TERMINADA COM FALHAS: {ok} de {n} genes concluídos, {failed} falharam",
        'en': "ANALYSIS FINISHED WITH FAILURES: {ok} of {n} genes completed, {failed} failed",
    },
    'summary_stopped': {
        'pt': "ANÁLISE INTERROMPIDA: {ok} de {n} genes concluídos antes de parar",
        'en': "ANALYSIS STOPPED: {ok} of {n} genes completed before stopping",
    },
    'summary_failed_item': {
        'pt': "  • {gene}: {reason}",
        'en': "  • {gene}: {reason}",
    },
    'results_in': {
        'pt': "Resultados em: {path}",
        'en': "Results in: {path}",
    },
}


def set_language(lang: str) -> None:
    global _LANG
    _LANG = 'pt' if str(lang).lower().startswith('pt') else 'en'


def get_language() -> str:
    return _LANG


def system_language() -> str:
    """'pt' se o sistema estiver em português, senão 'en'."""
    candidates = [os.environ.get(k, '') for k in ('LC_ALL', 'LC_MESSAGES', 'LANG', 'LANGUAGE')]
    try:
        loc = locale.getlocale()[0] or ''
        candidates.append(loc)
    except Exception:
        pass
    try:  # Windows: idioma da interface do usuário
        import ctypes
        lcid = ctypes.windll.kernel32.GetUserDefaultUILanguage()  # type: ignore[attr-defined]
        if (lcid & 0x3FF) == 0x16:  # LANG_PORTUGUESE
            return 'pt'
    except Exception:
        pass
    for c in candidates:
        if c and c.lower().startswith(('pt', 'portuguese')):
            return 'pt'
    return 'en'


def t(key: str, **kw) -> str:
    entry = MESSAGES.get(key)
    if not entry:
        return key
    template = entry.get(_LANG) or entry['en']
    try:
        return template.format(**kw)
    except (KeyError, IndexError, ValueError):
        return template
