"""
Peças de interface compartilhadas pela janela principal e pelo painel de
resultados: diálogos (validação antes de rodar, confirmação, Sobre), ajuste
de janela à tela, abrir pasta no gerenciador de arquivos, cores.
"""

import os
import platform
import subprocess
import sys
from pathlib import Path

import customtkinter as ctk

from .gui_texts import TEXTS, get_language

_ON_LINUX = platform.system() == "Linux"
_ON_WIN = platform.system() == "Windows"
FONT_UI = "DejaVu Sans" if _ON_LINUX else "Roboto"
FONT_MONO = "Cascadia Code" if _ON_WIN else "DejaVu Sans Mono"

# Paleta com contraste >= 4,5:1 (WCAG AA) para texto sobre os fundos escuros
# usados nos cartões (#0c0c0e .. #20202c). Ver tests/test_gui_contrast.py.
PALETTE = {
    'bg_dark': '#0d0d11',
    'bg_card': '#17171f',
    'bg_card_hover': '#20202c',
    'text_primary': '#eeeef2',
    'text_secondary': '#a3a3b8',
    'text_tertiary': '#8e8ea4',
    'text_muted': '#8a8aa0',
    'accent_text': '#818cf8',       # índigo legível como texto
    'accent_fill': '#4f46e5',       # índigo como fundo (texto branco 6,3:1)
    'success_text': '#22c55e',
    'success_fill': '#15803d',
    'warning_text': '#f59e0b',
    'warning_fill': '#b45309',
    'danger_text': '#f87171',
    'danger_fill': '#b91c1c',
    'info_text': '#22d3ee',
    'info_fill': '#0e7490',
    'border': '#262632',
    'border_hover': '#3a3a4e',
    # Camadas da janela principal (escuro -> claro), separadas por tom, não
    # por bordas. Texto primário/secundário/terciário >= 4,5:1 em todas.
    'bg_window': '#0a0b0e',          # fundo da janela
    'bg_panel': '#111318',           # painel lateral
    'bg_surface': '#181a21',         # cartões, abas, log
    'bg_elevated': '#20232c',        # cartão de modelo, botão de tom
    'bg_elevated_hover': '#292c37',  # hover sobre bg_elevated (só texto primário)
    'bg_inset': '#0e1014',           # campos de entrada
    'divider': '#262a34',            # linha de 1 px entre grupos
    'control_border': '#2c303c',     # contorno de botão secundário / campo
    'control_border_hover': '#3b4050',
    'neutral_fill': '#343846',       # botão neutro preenchido (texto branco 11,7:1)
    'switch_knob': '#e6e7ee',
    # Cores semânticas "subtle" (estilo Primer) do painel de resultados: fundo
    # quase neutro + texto claro só no valor/rótulo do veredito. Cada *_fg tem
    # >= 4,5:1 sobre o seu *_subtle e sobre todas as camadas acima.
    'success_subtle': '#0f2a1d',     # q < 0,05: fundo do rótulo "significativo"
    'success_fg': '#86efac',         # texto sobre success_subtle (10,9:1)
    'warning_subtle': '#2e2108',     # avisos (M7×M8 sem M8a×M8, análises órfãs)
    'warning_fg': '#fcd34d',         # texto sobre warning_subtle (10,9:1)
    'danger_subtle': '#341719',      # só falha de execução
    'danger_fg': '#fca5a5',          # texto sobre danger_subtle (8,6:1)
    'row_alt': '#15171d',            # zebra suave das tabelas (entre bg_panel e bg_surface)
}

# Escala de espaçamento (grade de 4 px) e raios de canto
SPACE = {'xs': 4, 'sm': 8, 'md': 12, 'lg': 16, 'xl': 24, 'xxl': 32}
RADIUS = {'field': 6, 'card': 8, 'panel': 12}
# Escala de fontes (pt): nada abaixo de 11
FONT_SIZE = {'xs': 11, 'sm': 12, 'md': 13, 'lg': 15, 'xl': 17, 'xxl': 20}


def mix(color_a: str, color_b: str, t: float) -> str:
    """Mistura duas cores hex: t = 0 -> a, t = 1 -> b."""
    a = [int(color_a.lstrip('#')[i:i + 2], 16) for i in (0, 2, 4)]
    b = [int(color_b.lstrip('#')[i:i + 2], 16) for i in (0, 2, 4)]
    return '#' + ''.join(f"{round(x + (y - x) * t):02x}" for x, y in zip(a, b))


def hover_tint(accent: str, background: str = '#17171f') -> str:
    """Fundo de hover de um botão de contorno: um tom escuro da cor de
    destaque, para o texto (na cor de destaque) continuar legível."""
    return mix(background, accent, 0.22)


def fit_to_screen(win, width: int, height: int, min_w: int = 1024, min_h: int = 640) -> None:
    """Abre a janela no tamanho pedido, mas nunca maior que a tela."""
    win.update_idletasks()
    sw, sh = win.winfo_screenwidth(), win.winfo_screenheight()
    w = min(width, sw - 40)
    h = min(height, sh - 80)
    x = max(0, (sw - w) // 2)
    y = max(0, (sh - h) // 3)
    win.geometry(f"{w}x{h}+{x}+{y}")
    win.minsize(min(min_w, w), min(min_h, h))


def disable_mouse_wheel(widget) -> None:
    """Impede que a roda do mouse mude o valor de um CTkSlider ao rolar o
    painel (a rolagem segue para o painel)."""
    canvas = getattr(widget, '_canvas', None)
    if canvas is None:
        return
    for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
        try:
            canvas.unbind(seq)
        except Exception:
            pass


def open_folder(path) -> bool:
    path = str(Path(path))
    try:
        if _ON_WIN:
            os.startfile(path)  # type: ignore[attr-defined]
        elif sys.platform == 'darwin':
            subprocess.Popen(['open', path])
        else:
            subprocess.Popen(['xdg-open', path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except Exception:
        return False


class _Modal(ctk.CTkToplevel):
    def __init__(self, parent, title: str, width: int, height: int):
        super().__init__(parent)
        self.title(title)
        self.configure(fg_color=PALETTE['bg_window'])
        fit_to_screen(self, width, height, min_w=420, min_h=240)
        self.transient(parent)
        self.bind("<Escape>", lambda e: self._close(None))
        self.protocol("WM_DELETE_WINDOW", lambda: self._close(None))
        self.result = None
        self.after(50, self._grab)

    def _grab(self):
        try:
            self.grab_set()
            self.focus_force()
        except Exception:
            pass

    def _close(self, value):
        self.result = value
        try:
            self.grab_release()
        except Exception:
            pass
        self.destroy()

    def show(self):
        self.wait_window()
        return self.result


def _button(parent, text, command, fill, **kw):
    return ctk.CTkButton(parent, text=text, command=command, fg_color=fill,
                         hover_color=mix(fill, '#000000', 0.2), text_color='#ffffff',
                         text_color_disabled=mix('#ffffff', fill, 0.45),
                         font=(FONT_UI, FONT_SIZE['md'], 'bold'), height=36,
                         corner_radius=RADIUS['card'], **kw)


def ask_yes_no(parent, title: str, message: str, yes: str = None, no: str = None) -> bool:
    """Confirmação com botões traduzidos (o messagebox do Tk usa Yes/No do sistema)."""
    dlg = _Modal(parent, title, 460, 220)
    ctk.CTkLabel(dlg, text=message, font=(FONT_UI, FONT_SIZE['md']), wraplength=400, justify='left',
                 text_color=PALETTE['text_primary']).pack(padx=SPACE['xl'], pady=(SPACE['xl'], SPACE['lg']),
                                                          anchor='w')
    row = ctk.CTkFrame(dlg, fg_color='transparent')
    row.pack(fill='x', padx=SPACE['xl'], pady=(0, SPACE['xl']), side='bottom')
    _button(row, no or TEXTS['btn_no'], lambda: dlg._close(False), PALETTE['neutral_fill']).pack(side='right')
    _button(row, yes or TEXTS['btn_yes'], lambda: dlg._close(True), PALETTE['accent_fill']).pack(
        side='right', padx=(0, SPACE['sm']))
    return bool(dlg.show())


def show_message(parent, title: str, message: str, kind: str = 'info') -> None:
    dlg = _Modal(parent, title, 560, 300)
    color = {'error': PALETTE['danger_text'], 'warning': PALETTE['warning_text']}.get(
        kind, PALETTE['text_primary'])
    # borda só em aviso/erro (comunica o tipo); mensagem comum separa pelo tom
    box = ctk.CTkTextbox(dlg, font=(FONT_UI, FONT_SIZE['md']), fg_color=PALETTE['bg_surface'],
                         text_color=color, wrap='word', corner_radius=RADIUS['card'],
                         border_width=1 if kind in ('error', 'warning') else 0,
                         border_color=mix(PALETTE['bg_surface'], color, 0.45))
    box.pack(fill='both', expand=True, padx=SPACE['xl'], pady=(SPACE['xl'], SPACE['md']))
    box.insert('end', message)
    box.configure(state='disabled')
    _button(dlg, "OK", lambda: dlg._close(True), PALETTE['accent_fill'], width=100).pack(
        pady=(0, SPACE['lg']))
    dlg.show()


class PreflightDialog(_Modal):
    """Problemas encontrados ANTES de rodar. Retorna 'continue' ou 'fix'."""

    def __init__(self, parent, report):
        super().__init__(parent, TEXTS['preflight_title'], 900, 620)
        lang = get_language()
        n_err = sum(1 for i in report.issues if i.severity == 'error')
        n_warn = sum(1 for i in report.issues if i.severity == 'warning')
        head = ctk.CTkFrame(self, fg_color='transparent')
        head.pack(fill='x', padx=SPACE['xl'], pady=(SPACE['lg'], SPACE['sm']))
        ctk.CTkLabel(head, text=TEXTS['preflight_heading'].format(
                         genes=len(report.genes), errors=n_err, warnings=n_warn),
                     font=(FONT_UI, FONT_SIZE['lg'], 'bold'), text_color=PALETTE['text_primary'],
                     justify='left', wraplength=840).pack(anchor='w')
        ctk.CTkLabel(head, text=TEXTS['preflight_explain'], font=(FONT_UI, FONT_SIZE['sm']),
                     text_color=PALETTE['text_secondary'], justify='left',
                     wraplength=840).pack(anchor='w', pady=(SPACE['xs'], 0))

        box = ctk.CTkTextbox(self, font=(FONT_MONO, FONT_SIZE['sm']), fg_color=PALETTE['bg_surface'],
                             text_color=PALETTE['text_primary'], wrap='word',
                             corner_radius=RADIUS['card'], border_width=1,
                             border_color=mix(PALETTE['bg_surface'], PALETTE['warning_text'], 0.35))
        box.pack(fill='both', expand=True, padx=SPACE['xl'], pady=SPACE['sm'])
        box.tag_config('error', foreground=PALETTE['danger_text'])
        box.tag_config('warning', foreground=PALETTE['warning_text'])
        box.tag_config('info', foreground=PALETTE['text_tertiary'])
        box.tag_config('gene', foreground=PALETTE['accent_text'])
        tags = {'error': TEXTS['preflight_tag_error'], 'warning': TEXTS['preflight_tag_warning'],
                'info': TEXTS['preflight_tag_info']}
        for gene, issues in report.by_gene().items():
            box.insert('end', (gene or TEXTS['preflight_general']) + "\n", 'gene')
            for issue in sorted(issues, key=lambda i: {'error': 0, 'warning': 1, 'info': 2}[i.severity]):
                box.insert('end', f"  [{tags[issue.severity]}] ", issue.severity)
                box.insert('end', issue.message(lang) + "\n")
            box.insert('end', "\n")
        box.configure(state='disabled')

        row = ctk.CTkFrame(self, fg_color='transparent')
        row.pack(fill='x', padx=SPACE['xl'], pady=(SPACE['sm'], SPACE['lg']))
        _button(row, TEXTS['preflight_btn_continue'], lambda: self._close('continue'),
                PALETTE['warning_fill']).pack(side='right')
        _button(row, TEXTS['preflight_btn_fix'], lambda: self._close('fix'),
                PALETTE['accent_fill']).pack(side='right', padx=(0, SPACE['sm']))
        ctk.CTkLabel(row, text=TEXTS['preflight_continue_hint'], font=(FONT_UI, FONT_SIZE['sm']),
                     text_color=PALETTE['text_tertiary'], wraplength=480,
                     justify='left').pack(side='left')


def show_about(parent) -> None:
    from backend.codeml_backend import codeml_version, find_codeml
    from backend.version import __version__
    path = find_codeml()
    ver = codeml_version(path) if path else None
    text = TEXTS['about_text'].format(
        version=__version__, codeml=path or TEXTS['about_codeml_missing'],
        codeml_version=ver or '?', python=platform.python_version(),
        platform=platform.platform())
    show_message(parent, TEXTS['about_title'], text)
