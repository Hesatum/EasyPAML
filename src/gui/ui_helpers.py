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
}


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
        self.configure(fg_color=PALETTE['bg_dark'])
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
                         font=(FONT_UI, 13, 'bold'), height=36, corner_radius=8, **kw)


def ask_yes_no(parent, title: str, message: str, yes: str = None, no: str = None) -> bool:
    """Confirmação com botões traduzidos (o messagebox do Tk usa Yes/No do sistema)."""
    dlg = _Modal(parent, title, 460, 220)
    ctk.CTkLabel(dlg, text=message, font=(FONT_UI, 13), wraplength=400, justify='left',
                 text_color=PALETTE['text_primary']).pack(padx=24, pady=(24, 16), anchor='w')
    row = ctk.CTkFrame(dlg, fg_color='transparent')
    row.pack(fill='x', padx=24, pady=(0, 20), side='bottom')
    _button(row, no or TEXTS['btn_no'], lambda: dlg._close(False), '#3d3d4a').pack(side='right')
    _button(row, yes or TEXTS['btn_yes'], lambda: dlg._close(True), PALETTE['accent_fill']).pack(
        side='right', padx=(0, 8))
    return bool(dlg.show())


def show_message(parent, title: str, message: str, kind: str = 'info') -> None:
    dlg = _Modal(parent, title, 560, 300)
    color = {'error': PALETTE['danger_text'], 'warning': PALETTE['warning_text']}.get(
        kind, PALETTE['text_primary'])
    box = ctk.CTkTextbox(dlg, font=(FONT_UI, 13), fg_color=PALETTE['bg_card'],
                         text_color=color, wrap='word', border_width=1,
                         border_color=PALETTE['border'])
    box.pack(fill='both', expand=True, padx=20, pady=(20, 10))
    box.insert('end', message)
    box.configure(state='disabled')
    _button(dlg, "OK", lambda: dlg._close(True), PALETTE['accent_fill'], width=100).pack(pady=(0, 16))
    dlg.show()


class PreflightDialog(_Modal):
    """Problemas encontrados ANTES de rodar. Retorna 'continue' ou 'fix'."""

    def __init__(self, parent, report):
        super().__init__(parent, TEXTS['preflight_title'], 900, 620)
        lang = get_language()
        n_err = sum(1 for i in report.issues if i.severity == 'error')
        n_warn = sum(1 for i in report.issues if i.severity == 'warning')
        head = ctk.CTkFrame(self, fg_color='transparent')
        head.pack(fill='x', padx=20, pady=(18, 6))
        ctk.CTkLabel(head, text=TEXTS['preflight_heading'].format(
                         genes=len(report.genes), errors=n_err, warnings=n_warn),
                     font=(FONT_UI, 15, 'bold'), text_color=PALETTE['text_primary'],
                     justify='left', wraplength=840).pack(anchor='w')
        ctk.CTkLabel(head, text=TEXTS['preflight_explain'], font=(FONT_UI, 12),
                     text_color=PALETTE['text_secondary'], justify='left',
                     wraplength=840).pack(anchor='w', pady=(4, 0))

        box = ctk.CTkTextbox(self, font=(FONT_MONO, 12), fg_color=PALETTE['bg_card'],
                             text_color=PALETTE['text_primary'], wrap='word',
                             border_width=1, border_color=PALETTE['border'])
        box.pack(fill='both', expand=True, padx=20, pady=8)
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
        row.pack(fill='x', padx=20, pady=(4, 18))
        _button(row, TEXTS['preflight_btn_continue'], lambda: self._close('continue'),
                PALETTE['warning_fill']).pack(side='right')
        _button(row, TEXTS['preflight_btn_fix'], lambda: self._close('fix'),
                PALETTE['accent_fill']).pack(side='right', padx=(0, 8))
        ctk.CTkLabel(row, text=TEXTS['preflight_continue_hint'], font=(FONT_UI, 12),
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
