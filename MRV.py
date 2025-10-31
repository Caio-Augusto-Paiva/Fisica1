import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import math

# utiliza o basico das teorias do MRV para resolver sistemas com 5 variaveis, onde 3 sao dadas e 2 sao desconhecidas
# a partir dai consegue resolver o sistema e mostrar o passo a passo


def solve_muv_system(s, v0, vf, a, t):

    def format_eq(equation_str, values_str, result_str):
        steps = ""
        steps += f"   Equação:   {equation_str}\n"
        steps += f"   Subst.:    {values_str}\n"
        steps += f"   Resultado: {result_str}\n\n"
        return steps

    steps = "--- Passo-a-Passo: MUV (Aceleração Constante) ---\n"

    try:
        vars = {'s': s, 'v0': v0, 'vf': vf, 'a': a, 't': t}

        unknown_keys = [k for k, v in vars.items() if v is None]
        known_vars = {k: v for k, v in vars.items() if v is not None}
        steps += f"Dados Fornecidos: {known_vars}\n"
        steps += f"A calcular: {unknown_keys}\n\n"

        for i in range(3):
            made_progress = False
            if all(v is not None for v in vars.values()):
                break

            # Eq 1: v = v0 + a*t
            if vars['vf'] is None and None not in [vars['v0'], vars['a'], vars['t']]:
                steps += f"Passo {i+1}: Resolvendo para 'vf' usando v = v0 + a*t\n"
                val = vars['v0'] + vars['a'] * vars['t']
                steps += format_eq("v = v0 + a * t",
                                   f"v = {vars['v0']} + {vars['a']} * {vars['t']}", f"v = {val:.2f}")
                vars['vf'] = val
                made_progress = True
            elif vars['v0'] is None and None not in [vars['vf'], vars['a'], vars['t']]:
                steps += f"Passo {i+1}: Resolvendo para 'v0' usando v = v0 + a*t\n"
                val = vars['vf'] - vars['a'] * vars['t']
                steps += format_eq("v0 = v - a * t",
                                   f"v0 = {vars['vf']} - {vars['a']} * {vars['t']}", f"v0 = {val:.2f}")
                vars['v0'] = val
                made_progress = True
            elif vars['a'] is None and None not in [vars['vf'], vars['v0'], vars['t']]:
                steps += f"Passo {i+1}: Resolvendo para 'a' usando v = v0 + a*t\n"
                if vars['t'] == 0:
                    raise ZeroDivisionError(
                        "Tempo 't' não pode ser 0 se 'a' é desconhecido.")
                val = (vars['vf'] - vars['v0']) / vars['t']
                steps += format_eq("a = (v - v0) / t",
                                   f"a = ({vars['vf']} - {vars['v0']}) / {vars['t']}", f"a = {val:.2f}")
                vars['a'] = val
                made_progress = True
            elif vars['t'] is None and None not in [vars['vf'], vars['v0'], vars['a']]:
                steps += f"Passo {i+1}: Resolvendo para 't' usando v = v0 + a*t\n"
                if vars['a'] == 0:
                    raise ZeroDivisionError(
                        "Aceleração 'a' não pode ser 0 se 't' é desconhecido (para esta eq.).")
                val = (vars['vf'] - vars['v0']) / vars['a']
                steps += format_eq("t = (v - v0) / a",
                                   f"t = ({vars['vf']} - {vars['v0']}) / {vars['a']}", f"t = {val:.2f}")
                vars['t'] = val
                made_progress = True

            # Eq 2: s = (v0 + v)t / 2
            elif vars['s'] is None and None not in [vars['v0'], vars['vf'], vars['t']]:
                steps += f"Passo {i+1}: Resolvendo para 's' usando s = (v0 + v)t / 2\n"
                val = (vars['v0'] + vars['vf']) / 2 * vars['t']
                steps += format_eq("s = [(v0 + v) / 2] * t",
                                   f"s = [({vars['v0']} + {vars['vf']}) / 2] * {vars['t']}", f"s = {val:.2f}")
                vars['s'] = val
                made_progress = True

            # Eq 3: v² = v0² + 2as (Torricelli)
            elif vars['vf'] is None and None not in [vars['v0'], vars['a'], vars['s']]:
                steps += f"Passo {i+1}: Resolvendo para 'vf' usando v² = v0² + 2as\n"
                val_sq = vars['v0']**2 + 2 * vars['a'] * vars['s']
                if val_sq < 0:
                    raise ValueError(
                        f"Valor negativo dentro da raiz (v0² + 2as = {val_sq:.2f}).")
                val = math.sqrt(val_sq)
                steps += format_eq("v² = v0² + 2 * a * s",
                                   f"v² = {vars['v0']}² + 2 * {vars['a']} * {vars['s']}", f"v² = {val_sq:.2f}  =>  v = {val:.2f}")
                vars['vf'] = val
                made_progress = True
            elif vars['v0'] is None and None not in [vars['vf'], vars['a'], vars['s']]:
                steps += f"Passo {i+1}: Resolvendo para 'v0' usando v² = v0² + 2as\n"
                val_sq = vars['vf']**2 - 2 * vars['a'] * vars['s']
                if val_sq < 0:
                    raise ValueError(
                        f"Valor negativo dentro da raiz (v² - 2as = {val_sq:.2f}).")
                val = math.sqrt(val_sq)
                steps += format_eq("v0² = v² - 2 * a * s",
                                   f"v0² = {vars['vf']}² - 2 * {vars['a']} * {vars['s']}", f"v0² = {val_sq:.2f}  =>  v0 = {val:.2f}")
                vars['v0'] = val
                made_progress = True
            elif vars['a'] is None and None not in [vars['vf'], vars['v0'], vars['s']]:
                steps += f"Passo {i+1}: Resolvendo para 'a' usando v² = v0² + 2as\n"
                if vars['s'] == 0:
                    raise ZeroDivisionError(
                        "Deslocamento 's' não pode ser 0 se 'a' é desconhecido.")
                val = (vars['vf']**2 - vars['v0']**2) / (2 * vars['s'])
                steps += format_eq("a = (v² - v0²) / (2s)",
                                   f"a = ({vars['vf']}² - {vars['v0']}²) / (2 * {vars['s']})", f"a = {val:.2f}")
                vars['a'] = val
                made_progress = True
            elif vars['s'] is None and None not in [vars['vf'], vars['v0'], vars['a']]:
                steps += f"Passo {i+1}: Resolvendo para 's' usando v² = v0² + 2as\n"
                if vars['a'] == 0:
                    raise ZeroDivisionError(
                        "Aceleração 'a' não pode ser 0 se 's' é desconhecido.")
                val = (vars['vf']**2 - vars['v0']**2) / (2 * vars['a'])
                steps += format_eq("s = (v² - v0²) / (2a)",
                                   f"s = ({vars['vf']}² - {vars['v0']}²) / (2 * {vars['a']})", f"s = {val:.2f}")
                vars['s'] = val
                made_progress = True

            # Eq 4: s = v0*t + 0.5*a*t²
            elif vars['s'] is None and None not in [vars['v0'], vars['a'], vars['t']]:
                steps += f"Passo {i+1}: Resolvendo para 's' usando s = v0*t + 0.5*a*t²\n"
                val = vars['v0'] * vars['t'] + 0.5 * vars['a'] * vars['t']**2
                steps += format_eq("s = v0*t + 0.5*a*t²",
                                   f"s = {vars['v0']}*{vars['t']} + 0.5*{vars['a']}*{vars['t']}²", f"s = {val:.2f}")
                vars['s'] = val
                made_progress = True

            if not made_progress and any(v is None for v in vars.values()):
                raise ValueError(
                    "Sistema insolúvel. Equações restantes dependem uma da outra.")

        if any(v is None for v in vars.values()):
            raise ValueError(
                f"Não foi possível resolver o sistema. Variáveis restantes: {[k for k, v in vars.items() if v is None]}")

        steps += "\n--- RESULTADO FINAL (Todos os valores) ---\n"
        steps += f"   Deslocamento (s):      {vars['s']:.2f} m\n"
        steps += f"   Velocidade Inicial (v0): {vars['v0']:.2f} m/s\n"
        steps += f"   Velocidade Final (vf):   {vars['vf']:.2f} m/s\n"
        steps += f"   Aceleração (a):        {vars['a']:.2f} m/s²\n"
        steps += f"   Tempo (t):             {vars['t']:.2f} s\n"

        return steps
    except Exception as e:
        return f"Erro no cálculo: {e}\n\nVerifique se os dados são fisicamente possíveis e se não há divisão por zero (ex: 't' ou 'a' = 0 quando são divisores)."


def solve_calculus_kinematics(coeffs_str, mode, const_str):
    try:
        if not coeffs_str:
            raise ValueError("Campo 'Coeficientes' está vazio.")
        coeffs_list = [float(c.strip()) for c in coeffs_str.split(',')]
        C1 = float(const_str or 0)
        poly_in = np.poly1d(coeffs_list)
    except Exception as e:
        raise ValueError(f"Erro ao ler coeficientes/constante: {e}")

    steps = "--- Passo-a-Passo: Cinemática (Cálculo com NumPy) ---\n"

    if mode == "Posição s(t)":
        s_poly = poly_in
        steps += f"1. Função Posição s(t) fornecida:\n(Usando np.poly1d: {coeffs_list})\n"
        steps += f"   s(t) =\n{s_poly}\n\n"
        v_poly = s_poly.deriv()
        steps += f"2. Derivando s(t) para encontrar v(t) [s.deriv()]:\n"
        steps += f"   v(t) = ds/dt =\n{v_poly}\n\n"
        a_poly = v_poly.deriv()
        steps += f"3. Derivando v(t) para encontrar a(t) [v.deriv()]:\n"
        steps += f"   a(t) = dv/dt =\n{a_poly}\n\n"

    elif mode == "Velocidade v(t)":
        v_poly = poly_in
        steps += f"1. Função Velocidade v(t) fornecida:\n(Usando np.poly1d: {coeffs_list})\n"
        steps += f"   v(t) =\n{v_poly}\n\n"
        s_poly = v_poly.integ(k=C1)
        steps += f"2. Integrando v(t) para encontrar s(t) [v.integ(k={C1})]:\n"
        steps += f"   (Constante C1 = s(0) = {C1})\n"
        steps += f"   s(t) = ∫v(t) dt =\n{s_poly}\n\n"
        a_poly = v_poly.deriv()
        steps += f"3. Derivando v(t) para encontrar a(t) [v.deriv()]:\n"
        steps += f"   a(t) = dv/dt =\n{a_poly}\n\n"

    elif mode == "Aceleração a(t)":
        a_poly = poly_in
        steps += f"1. Função Aceleração a(t) fornecida:\n(Usando np.poly1d: {coeffs_list})\n"
        steps += f"   a(t) =\n{a_poly}\n\n"
        v_poly = a_poly.integ(k=C1)
        steps += f"2. Integrando a(t) para encontrar v(t) [a.integ(k={C1})]:\n"
        steps += f"   (Constante C1 = v(0) = {C1})\n"
        steps += f"   v(t) = ∫a(t) dt =\n{v_poly}\n\n"
        s_poly = v_poly.integ(k=0)
        steps += f"3. Integrando v(t) para encontrar s(t) [v.integ(k=0)]:\n"
        steps += f"   (Assumindo Constante C2 = s(0) = 0)\n"
        steps += f"   s(t) = ∫v(t) dt =\n{s_poly}\n\n"

    steps += "--- RESULTADO FINAL ---\n"
    if 's_poly' in locals():
        steps += f"s(t) =\n{s_poly}\n\n"
    if 'v_poly' in locals():
        steps += f"v(t) =\n{v_poly}\n\n"
    if 'a_poly' in locals():
        steps += f"a(t) =\n{a_poly}\n"

    return steps

# --- Classe da Aplicação GUI ---


class KinematicsCalculatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Calculadora de Cinemática Unificada")
        self.root.geometry("700x550")

        # Dicionários para guardar widgets
        self.muv_entries = {}
        self.calc_widgets = {}

        # Variável para controlar o modo
        self.calc_mode = tk.StringVar(value="MUV")

        # --- Layout Principal ---
        main_frame = ttk.Frame(root, padding="10")
        main_frame.pack(fill='both', expand=True)

        # 1. Frame de Seleção de Modo
        self.create_mode_selector(main_frame)

        # 2. Frame de Entradas (que mudará)
        # Criamos um "container" para os frames de entrada
        self.input_container = ttk.Frame(main_frame)
        self.input_container.pack(side='left', fill='y', padx=5, pady=5)

        self.muv_input_frame = ttk.LabelFrame(
            self.input_container, text="Entradas MUV", padding="10")
        self.calculus_input_frame = ttk.LabelFrame(
            self.input_container, text="Entradas Cálculo", padding="10")

        self.create_muv_inputs(self.muv_input_frame)
        self.create_calculus_inputs(self.calculus_input_frame)

        # 3. Frame de Saída (Único)
        output_frame = ttk.LabelFrame(
            main_frame, text="Passo a Passo e Resultados", padding="10")
        output_frame.pack(side='right', fill='both',
                          expand=True, padx=5, pady=5)

        scrollbar = ttk.Scrollbar(output_frame, orient='vertical')
        self.output_text = tk.Text(output_frame, wrap='word', height=15, font=(
            "Courier New", 10), yscrollcommand=scrollbar.set)
        scrollbar.config(command=self.output_text.yview)
        scrollbar.pack(side='right', fill='y')
        self.output_text.pack(fill='both', expand=True)

        # 4. Botão de Cálculo (Único)
        calculate_button = ttk.Button(
            self.input_container, text="Calcular", command=self.perform_calculation)
        calculate_button.pack(side='bottom', fill='x', pady=10, ipady=5)

        # Configuração inicial
        self.update_interface()  # Mostra o frame correto
        self.set_result_text("Selecione o modo e preencha os dados.")

    def create_mode_selector(self, parent):
        frame = ttk.LabelFrame(parent, text="Modo de Cálculo", padding="10")
        frame.pack(fill='x', pady=5)

        rb1 = ttk.Radiobutton(frame, text="Aceleração Constante (MUV)",
                              variable=self.calc_mode, value="MUV", command=self.update_interface)
        rb1.pack(side='left', padx=10)

        rb2 = ttk.Radiobutton(frame, text="Aceleração Variável (Cálculo)",
                              variable=self.calc_mode, value="Calculus", command=self.update_interface)
        rb2.pack(side='left', padx=10)

    def create_muv_inputs(self, parent_frame):
        """Cria os widgets de entrada para MUV."""
        ttk.Label(parent_frame, text="Preencha EXATAMENTE 3 dos 5 campos.",
                  wraplength=180).pack(pady=(0, 10))
        self.muv_entries['s'] = self.create_entry_pair(
            parent_frame, "Deslocamento (s) [m]:")
        self.muv_entries['v0'] = self.create_entry_pair(
            parent_frame, "Velocidade Inicial (v0) [m/s]:")
        self.muv_entries['vf'] = self.create_entry_pair(
            parent_frame, "Velocidade Final (vf) [m/s]:")
        self.muv_entries['a'] = self.create_entry_pair(
            parent_frame, "Aceleração (a) [m/s²]:")
        self.muv_entries['t'] = self.create_entry_pair(
            parent_frame, "Tempo (t) [s]:")

    def create_calculus_inputs(self, parent_frame):
        """Cria os widgets de entrada para Cálculo."""
        ttk.Label(parent_frame, text="Função de Entrada:",
                  anchor='w').pack(fill='x')
        self.calc_widgets['mode'] = ttk.Combobox(parent_frame, values=[
                                                 "Posição s(t)", "Velocidade v(t)", "Aceleração a(t)"], state='readonly', width=25)
        self.calc_widgets['mode'].pack(fill='x', padx=5, pady=5)
        self.calc_widgets['mode'].current(0)

        self.calc_widgets['coeffs'] = self.create_entry_pair(
            parent_frame, "Coeficientes (ex: 5, 0, 2):")
        ttk.Label(parent_frame, text="(Para 5t² + 0t + 2)",
                  font=('Arial', 8)).pack(anchor='w', padx=5)

        self.calc_widgets['const'] = self.create_entry_pair(
            parent_frame, "Constante Integração (C):", pady=(10, 2))
        ttk.Label(parent_frame, text="(Para s(0) ou v(0). Deixe vazio para 0)", font=(
            'Arial', 8)).pack(anchor='w', padx=5)

    def update_interface(self):
        """Mostra/oculta os frames de entrada com base no modo selecionado."""
        mode = self.calc_mode.get()

        if mode == "MUV":
            # Oculta o frame de cálculo
            self.calculus_input_frame.pack_forget()
            # Mostra o frame MUV
            self.muv_input_frame.pack(fill='x', expand=True)
            self.set_result_text(
                "Deixe 2 campos de MUV vazios para serem calculados.")

        elif mode == "Calculus":
            # Oculta o frame MUV
            self.muv_input_frame.pack_forget()
            # Mostra o frame de cálculo
            self.calculus_input_frame.pack(fill='x', expand=True)
            self.set_result_text(
                "Insira os coeficientes do polinômio (do maior grau para o menor) separados por vírgula.")

    # --- Funções Auxiliares da GUI ---

    def create_entry_pair(self, parent, label_text, pady=2):
        frame = ttk.Frame(parent)
        label = ttk.Label(frame, text=label_text, width=25)
        label.pack(side='left', padx=5, pady=pady)
        entry = ttk.Entry(frame, width=15)
        entry.pack(side='left', padx=5, pady=pady)
        frame.pack(anchor='w')
        return entry

    def set_result_text(self, text):
        self.output_text.config(state='normal')
        self.output_text.delete('1.0', 'end')
        self.output_text.insert('1.0', text)
        self.output_text.config(state='disabled')

    def get_float_or_none(self, entry):
        val_str = entry.get().replace(',', '.')
        if not val_str:
            return None
        try:
            return float(val_str)
        except ValueError:
            raise ValueError(f"Entrada inválida: '{val_str}' não é um número.")

    def perform_calculation(self):
        """Callback único para o botão 'Calcular'."""
        mode = self.calc_mode.get()

        try:
            if mode == "MUV":
                s = self.get_float_or_none(self.muv_entries['s'])
                v0 = self.get_float_or_none(self.muv_entries['v0'])
                vf = self.get_float_or_none(self.muv_entries['vf'])
                a = self.get_float_or_none(self.muv_entries['a'])
                t = self.get_float_or_none(self.muv_entries['t'])

                inputs = [s, v0, vf, a, t]
                none_count = sum(1 for x in inputs if x is None)

                if none_count != 2:
                    raise ValueError(
                        f"Você deve preencher exatamente 3 campos (foram preenchidos {5-none_count}).")

                steps = solve_muv_system(s, v0, vf, a, t)
                self.set_result_text(steps)

            elif mode == "Calculus":
                mode_func = self.calc_widgets['mode'].get()
                coeffs_str = self.calc_widgets['coeffs'].get()
                const_str = self.calc_widgets['const'].get()

                steps = solve_calculus_kinematics(
                    coeffs_str, mode_func, const_str)
                self.set_result_text(steps)

        except Exception as e:
            messagebox.showerror("Erro de Entrada ou Cálculo", str(e))
            self.set_result_text(f"Erro: {e}")


if __name__ == "__main__":
    main_window = tk.Tk()
    app = KinematicsCalculatorApp(main_window)
    main_window.mainloop()
