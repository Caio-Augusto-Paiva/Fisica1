import tkinter as tk
from tkinter import ttk, messagebox
import math
import numpy as np  # <-- Importado NumPy

# --- Função de Lógica de Cálculo (com NumPy) ---


def calculate_2d_collision_unified(m1, v1ix, v1iy, m2, v2ix, v2iy, e):
    """
    Calcula uma colisão 2D usando vetores NumPy.
    - Se e = 0, trata como perfeitamente inelástica (grudam em 2D).
    - Se e > 0, assume que o impacto (impulso) ocorre ao longo do eixo X.
    """

    if e < 0 or e > 1:
        raise ValueError("Coeficiente 'e' deve estar entre 0 e 1.")

    # --- 1. Análise Inicial (com Vetores) ---
    steps = "--- 1. Análise Inicial do Sistema (Vetorial) ---\n"

    # Criando vetores de velocidade inicial
    v1i_vec = np.array([v1ix, v1iy])
    v2i_vec = np.array([v2ix, v2iy])
    steps += "Representando velocidades como vetores [vx, vy]:\n"
    steps += f"   v1i = {v1i_vec} m/s\n"
    steps += f"   v2i = {v2i_vec} m/s\n"

    # Momento Vetorial Inicial (Pi_vec)
    Pi_vec = (m1 * v1i_vec) + (m2 * v2i_vec)
    steps += "\nMomento Vetorial Inicial (Pi = m1*v1i + m2*v2i):\n"
    steps += f"   Pi = {m1} * {v1i_vec} + {m2} * {v2i_vec}\n"
    steps += f"   Pi = {Pi_vec}  (kg·m/s)  [Pix, Piy]\n"

    # Energia Cinética Inicial (Ki)
    # Ki = 0.5 * m * |v|²  ou  0.5 * m * (v · v)
    v1i_mag_sq = np.dot(v1i_vec, v1i_vec)  # np.dot(v,v) é |v|²
    v2i_mag_sq = np.dot(v2i_vec, v2i_vec)
    Ki = 0.5 * m1 * v1i_mag_sq + 0.5 * m2 * v2i_mag_sq
    steps += "\nEnergia Cinética Inicial (Ki = 0.5*m*(v·v)):\n"
    steps += f"   Ki = 0.5*{m1}*({v1i_mag_sq:.2f}) + 0.5*{m2}*({v2i_mag_sq:.2f})\n"
    steps += f"   Ki = {Ki:.2f} J\n"

    M_total = m1 + m2
    if M_total <= 0:
        raise ValueError("A soma das massas deve ser positiva.")

    # Inicializa vetores finais
    v1f_vec = np.zeros(2)
    v2f_vec = np.zeros(2)

    # --- LÓGICA DE CÁLCULO ---

    if e == 0:
        # CASO 1: Colisão Perfeitamente Inelástica (e=0)
        # Lógica vetorial é muito limpa aqui

        steps += f"\n--- 2. Cálculo (Colisão Perfeitamente Inelástica, e=0) ---\n"
        steps += "Os corpos grudam. O momento vetorial se conserva.\n"
        steps += "   Pi_vec = (m1 + m2) * vf_vec\n"

        vf_vec = Pi_vec / M_total  # vf_vec é um vetor [vfx, vfy]

        steps += f"   vf_vec = Pi_vec / M_total = {Pi_vec} / {M_total}\n"
        steps += f"   vf_vec = {vf_vec} m/s\n"

        v1f_vec = vf_vec
        v2f_vec = vf_vec

    else:
        # CASO 2: Colisão Parcial ou Elástica (e > 0)
        # Suposição: O impacto (força/impulso) ocorre SOMENTE ao longo do eixo X.

        if e == 1:
            steps += f"\n--- 2. Cálculo (Colisão Perfeitamente Elástica, e=1) ---\n"
        else:
            steps += f"\n--- 2. Cálculo (Colisão Parcialmente Elástica, e={e}) ---\n"

        steps += "[Suposição: Impacto (impulso) ocorre ao longo do Eixo X]\n"

        # Cálculo em Y (Componentes Y dos vetores não mudam)
        steps += f"Componentes Y (tangenciais) não são alteradas:\n"
        v1fy = v1iy
        v2fy = v2iy
        steps += f"   v1fy = v1iy = {v1fy:.2f} m/s\n"
        steps += f"   v2fy = v2iy = {v2fy:.2f} m/s\n"

        # Cálculo em X (Tratado como uma colisão 1D nos componentes X)
        steps += f"\nResolvendo para os Componentes X (Conservação de Px e 'e'):\n"
        Pix = Pi_vec[0]  # Pega o componente X do vetor Pi
        steps += f"   (Eq 1) Px: {Pix:.2f} = {m1}*v1fx + {m2}*v2fx\n"
        steps += f"   (Eq 2) 'e': v2fx - v1fx = e * (v1ix - v2ix)\n"

        v_rel_ix = v1ix - v2ix  # Velocidade relativa inicial em X

        v1fx = (Pix - m2 * e * v_rel_ix) / M_total
        v2fx = v1fx + e * v_rel_ix

        steps += f"Resolvendo o sistema para os componentes X:\n"
        steps += f"   v1fx = {v1fx:.2f} m/s\n"
        steps += f"   v2fx = {v2fx:.2f} m/s\n"

        # Monta os vetores finais
        v1f_vec = np.array([v1fx, v1fy])
        v2f_vec = np.array([v2fx, v2fy])

    # --- 3. ANÁLISE FINAL (Vetorial) ---

    steps += "\n--- 3. Análise Final do Sistema ---\n"
    steps += "Resultados das Velocidades Finais (Vetores):\n"
    # np.round(..., 2) arredonda os elementos do array
    steps += f"   v1f_vec = {np.round(v1f_vec, 2)} m/s\n"
    steps += f"   v2f_vec = {np.round(v2f_vec, 2)} m/s\n"

    # Energia Cinética Final (Kf)
    v1f_mag_sq = np.dot(v1f_vec, v1f_vec)
    v2f_mag_sq = np.dot(v2f_vec, v2f_vec)
    Kf = (0.5 * m1 * v1f_mag_sq) + (0.5 * m2 * v2f_mag_sq)
    steps += "\nEnergia Cinética Final (Kf = 0.5*m*(vf·vf)):\n"
    steps += f"   Kf = 0.5*{m1}*({v1f_mag_sq:.2f}) + 0.5*{m2}*({v2f_mag_sq:.2f})\n"
    steps += f"   Kf = {Kf:.2f} J\n"

    # Perda de Energia (Delta K)
    delta_K = Kf - Ki
    steps += f"\nVariação de Energia (ΔK = Kf - Ki): {delta_K:.2f} J\n"

    if abs(delta_K) < 1e-6:  # Comparação segura para floats
        steps += "   Tipo: Perfeitamente Elástica (Energia conservada).\n"
    elif e == 0:
        steps += "   Tipo: Perfeitamente Inelástica (Máxima perda de energia cinética).\n"
    else:
        steps += f"   Tipo: Parcialmente Elástica (Energia perdida).\n"

    return steps

# --- Classe da Aplicação GUI ---
# (O código da GUI não precisa de NENHUMA alteração)
# Ele continua coletando os mesmos 7 escalares (m1, v1ix, v1iy, m2, v2ix, v2iy, e)
# e passando para a função de lógica.


class CollisionCalculatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Calculadora de Colisão Unificada (2D c/ NumPy)")
        self.root.geometry("800x600")

        self.entries = {}

        container = ttk.Frame(root)
        container.pack(fill='both', expand=True, padx=10, pady=10)

        input_frame = ttk.LabelFrame(
            container, text="Dados de Entrada", padding="10")
        input_frame.pack(side='left', fill='y', padx=5, pady=5)

        output_frame = ttk.LabelFrame(
            container, text="Passo a Passo e Resultados", padding="10")
        output_frame.pack(side='right', fill='both',
                          expand=True, padx=5, pady=5)

        self.create_input_frame(input_frame)
        self.create_output_frame(output_frame)

    def create_entry_pair(self, parent, entry_key, label_text):
        frame = ttk.Frame(parent)
        label = ttk.Label(frame, text=label_text, width=25)
        label.pack(side='left', padx=5, pady=3)
        entry = ttk.Entry(frame, width=15)
        entry.pack(side='left', padx=5, pady=3)
        frame.pack(anchor='w')
        self.entries[entry_key] = entry

    def create_input_frame(self, parent_frame):
        ttk.Label(parent_frame, text="Corpo 1:", font=(
            'Arial', 10, 'bold')).pack(anchor='w')
        self.create_entry_pair(parent_frame, 'm1', "Massa (m1) [kg]:")
        self.create_entry_pair(parent_frame, 'v1ix',
                               "Velocidade Inicial X (v1ix):")
        self.create_entry_pair(parent_frame, 'v1iy',
                               "Velocidade Inicial Y (v1iy):")

        ttk.Label(parent_frame, text="Corpo 2:", font=(
            'Arial', 10, 'bold')).pack(anchor='w', pady=(15, 0))
        self.create_entry_pair(parent_frame, 'm2', "Massa (m2) [kg]:")
        self.create_entry_pair(parent_frame, 'v2ix',
                               "Velocidade Inicial X (v2ix):")
        self.create_entry_pair(parent_frame, 'v2iy',
                               "Velocidade Inicial Y (v2iy):")

        ttk.Label(parent_frame, text="Dados da Colisão:", font=(
            'Arial', 10, 'bold')).pack(anchor='w', pady=(15, 0))
        self.create_entry_pair(
            parent_frame, 'e', "Coeficiente 'e':")

        calculate_button = ttk.Button(
            parent_frame, text="Calcular Colisão", command=self.perform_calculation)
        calculate_button.pack(side='bottom', pady=20, fill='x', ipady=5)

    def create_output_frame(self, parent_frame):
        scrollbar = ttk.Scrollbar(parent_frame, orient='vertical')
        self.result_text = tk.Text(parent_frame, wrap='word', height=15, font=(
            "Courier New", 10), yscrollcommand=scrollbar.set)
        scrollbar.config(command=self.result_text.yview)

        scrollbar.pack(side='right', fill='y')
        self.result_text.pack(fill='both', expand=True)

        self.set_result_text("Preencha todos os dados e clique em 'Calcular Colisão'.\n\n"
                             "Lembretes:\n"
                             "- Campos de velocidade vazios serão tratados como 0.\n"
                             "- 'm1', 'm2' e 'e' são obrigatórios.\n"
                             "- e=0: Colisão perfeitamente inelástica (grudam).\n"
                             "- e=1: Colisão perfeitamente elástica.")

    def set_result_text(self, text):
        self.result_text.config(state='normal')
        self.result_text.delete('1.0', 'end')
        self.result_text.insert('1.0', text)
        self.result_text.config(state='disabled')

    def get_float(self, entry_key):
        is_velocity_field = entry_key.startswith('v')

        try:
            value_str = self.entries[entry_key].get().replace(',', '.')

            if not value_str:
                if is_velocity_field:
                    return 0.0
                else:
                    raise ValueError(
                        f"Campo '{entry_key}' é obrigatório e está vazio.")

            return float(value_str)
        except ValueError as ve:
            if "obrigatório" in str(ve):
                raise ve
            raise ValueError(
                f"Entrada inválida para '{entry_key}'. Use apenas números.")

    def perform_calculation(self):
        try:
            # Coleta de dados (continua igual)
            m1 = self.get_float('m1')
            v1ix = self.get_float('v1ix')
            v1iy = self.get_float('v1iy')
            m2 = self.get_float('m2')
            v2ix = self.get_float('v2ix')
            v2iy = self.get_float('v2iy')
            e = self.get_float('e')

            # Calcular (chama a nova função baseada em NumPy)
            steps = calculate_2d_collision_unified(
                m1, v1ix, v1iy, m2, v2ix, v2iy, e)

            self.set_result_text(steps)

        except ValueError as ve:
            messagebox.showerror("Erro de Entrada", str(ve))
            self.set_result_text(f"Erro: {ve}")
        except Exception as e:
            messagebox.showerror(
                "Erro de Cálculo", f"Ocorreu um erro inesperado: {e}")
            self.set_result_text(f"Erro: {e}")


# --- Ponto de Entrada Principal ---
if __name__ == "__main__":
    main_window = tk.Tk()
    app = CollisionCalculatorApp(main_window)
    main_window.mainloop()
