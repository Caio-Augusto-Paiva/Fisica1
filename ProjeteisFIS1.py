import tkinter as tk
from tkinter import ttk, messagebox
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# --- Constantes ---
GRAVITY = 9.81  # Aceleração da gravidade (m/s^2)

# --- Funções de Lógica de Cálculo (AGORA RETORNAM 'steps') ---


def calculate_range_and_time(v0, theta_deg, y0, y_final):
    """
    Calcula o alcance (x) e o tempo total de voo (t).
    Retorna: (range_x, time, steps_string)
    """
    steps = f"--- Passo-a-passo: Alcance e Tempo de Voo ---\n"
    steps += f"1. Dados Fornecidos:\n"
    steps += f"   Velocidade Inicial (v0): {v0} m/s\n"
    steps += f"   Ângulo (θ): {theta_deg}°\n"
    steps += f"   Altura Inicial (y0): {y0} m\n"
    steps += f"   Altura Final (y_final): {y_final} m\n"
    steps += f"   Gravidade (g): {GRAVITY} m/s²\n"

    theta_rad = math.radians(theta_deg)
    v0y = v0 * math.sin(theta_rad)
    v0x = v0 * math.cos(theta_rad)
    steps += f"\n2. Decomposição da Velocidade Inicial:\n"
    steps += f"   v0x = v0 * cos(θ) = {v0} * cos({theta_deg}°) = {v0x:.2f} m/s\n"
    steps += f"   v0y = v0 * sin(θ) = {v0} * sin({theta_deg}°) = {v0y:.2f} m/s\n"

    steps += f"\n3. Encontrando o Tempo de Voo (t):\n"
    steps += f"   Usamos a Equação da Posição Vertical:\n"
    steps += f"   y(t) = y0 + v0y*t - 0.5*g*t²\n"
    steps += f"   Substituindo: {y_final} = {y0} + ({v0y:.2f})*t - 0.5*({GRAVITY})*t²\n"
    steps += f"   Reorganizando (forma a*t² + b*t + c = 0):\n"

    a = 0.5 * GRAVITY
    b = -v0y
    c = y_final - y0
    steps += f"   ({a:.2f})t² + ({b:.2f})t + ({c:.2f}) = 0\n"

    discriminant = b**2 - 4*a*c
    steps += f"   Resolvendo com Bhaskara (Δ = b² - 4ac):\n"
    steps += f"   Δ = ({b:.2f})² - 4*({a:.2f})*({c:.2f}) = {discriminant:.2f}\n"

    if discriminant < 0:
        raise ValueError("Trajetória impossível (discriminante negativo).")

    t1 = (-b + math.sqrt(discriminant)) / (2*a)
    t2 = (-b - math.sqrt(discriminant)) / (2*a)
    time = max(t1, t2)
    steps += f"   t = (-b ± √Δ) / 2a  => t1 = {t1:.2f} s, t2 = {t2:.2f} s\n"
    steps += f"   O tempo de voo é a maior raiz positiva: t = {time:.2f} s\n"

    if time < 0:
        raise ValueError("Tempo de voo negativo. Verifique as entradas.")

    steps += f"\n4. Encontrando o Alcance (x):\n"
    steps += f"   Usamos a Equação da Posição Horizontal:\n"
    steps += f"   x = v0x * t\n"
    range_x = v0x * time
    steps += f"   x = {v0x:.2f} * {time:.2f} = {range_x:.2f} m\n"

    steps += f"\n--- RESULTADO FINAL ---\n"
    steps += f"   Alcance Total (x): {range_x:.2f} m\n"
    steps += f"   Tempo de Voo (t): {time:.2f} s\n"

    return range_x, time, steps


def calculate_max_height(v0, theta_deg, y0):
    """
    Calcula a altura máxima (ymax) e o tempo para atingi-la (t_up).
    Retorna: (max_height, time_to_peak, steps_string)
    """
    steps = f"--- Passo-a-passo: Altura Máxima ---\n"
    steps += f"1. Dados Fornecidos:\n"
    steps += f"   Velocidade Inicial (v0): {v0} m/s\n"
    steps += f"   Ângulo (θ): {theta_deg}°\n"
    steps += f"   Altura Inicial (y0): {y0} m\n"

    theta_rad = math.radians(theta_deg)
    v0y = v0 * math.sin(theta_rad)
    steps += f"\n2. Decomposição da Velocidade Vertical:\n"
    steps += f"   v0y = v0 * sin(θ) = {v0} * sin({theta_deg}°) = {v0y:.2f} m/s\n"

    steps += f"\n3. Encontrando a Altura Máxima (y_max):\n"
    steps += f"   No ponto mais alto, a velocidade vertical (vy) é 0.\n"
    steps += f"   Usamos a Equação de Torricelli:\n"
    steps += f"   vy² = v0y² - 2*g*(y_max - y0)\n"
    steps += f"   0² = ({v0y:.2f})² - 2*({GRAVITY})*(y_max - {y0})\n"
    steps += f"   Isolando y_max:\n"
    steps += f"   y_max = y0 + (v0y² / (2*g))\n"
    max_height = y0 + (v0y**2) / (2 * GRAVITY)
    steps += f"   y_max = {y0} + (({v0y:.2f})² / (2 * {GRAVITY})) = {max_height:.2f} m\n"

    steps += f"\n4. Encontrando o Tempo para Atingir o Pico (t_pico):\n"
    steps += f"   Usamos a Equação da Velocidade Vertical (vy = 0):\n"
    steps += f"   vy = v0y - g*t\n"
    steps += f"   0 = {v0y:.2f} - {GRAVITY} * t_pico\n"
    time_to_peak = v0y / GRAVITY
    steps += f"   t_pico = {v0y:.2f} / {GRAVITY} = {time_to_peak:.2f} s\n"

    steps += f"\n--- RESULTADO FINAL ---\n"
    steps += f"   Altura Máxima (y_max): {max_height:.2f} m\n"
    steps += f"   Tempo até o Pico (t_pico): {time_to_peak:.2f} s\n"

    return max_height, time_to_peak, steps


def calculate_initial_velocity(theta_deg, y0, y_final, range_x):
    """
    Calcula a velocidade inicial (v0) necessária.
    Retorna: (v0, steps_string)
    """
    steps = f"--- Passo-a-passo: Velocidade Inicial ---\n"
    steps += f"1. Dados Fornecidos:\n"
    steps += f"   Ângulo (θ): {theta_deg}°\n"
    steps += f"   Altura Inicial (y0): {y0} m\n"
    steps += f"   Altura Final (y_final): {y_final} m\n"
    steps += f"   Alcance (x): {range_x} m\n"

    theta_rad = math.radians(theta_deg)
    tan_theta = math.tan(theta_rad)
    cos_theta = math.cos(theta_rad)

    steps += f"\n2. Usando a Equação da Trajetória (sem tempo):\n"
    steps += f"   y = y0 + x*tan(θ) - (g*x²) / (2*v0²*cos²(θ))\n"
    steps += f"   Isolando v0²:\n"
    steps += f"   (g*x²) / (2*v0²*cos²(θ)) = y0 - y_final + x*tan(θ)\n"
    steps += f"   v0² = (g*x²) / [ 2*cos²(θ) * (y0 - y_final + x*tan(θ)) ]\n"

    steps += f"\n3. Substituindo os valores:\n"
    numerator = GRAVITY * range_x**2
    denominator_part = y0 - y_final + range_x * tan_theta
    denominator = 2 * (cos_theta**2) * denominator_part

    steps += f"   Numerador = {GRAVITY} * {range_x}² = {numerator:.2f}\n"
    steps += f"   Denominador = 2 * cos²({theta_deg}) * ({y0} - {y_final} + {range_x}*tan({theta_deg}))\n"
    steps += f"   Denominador = 2 * ({cos_theta:.2f})² * ({denominator_part:.2f}) = {denominator:.2f}\n"

    if denominator <= 0:
        raise ValueError(
            "Trajetória impossível com este ângulo (alvo inatingível).")

    v0_squared = numerator / denominator
    v0 = math.sqrt(v0_squared)
    steps += f"   v0² = {numerator:.2f} / {denominator:.2f} = {v0_squared:.2f}\n"
    steps += f"   v0 = √{v0_squared:.2f} = {v0:.2f} m/s\n"

    steps += f"\n--- RESULTADO FINAL ---\n"
    steps += f"   Velocidade Inicial (v0): {v0:.2f} m/s\n"

    return v0, steps


def calculate_angle(v0, y0, y_final, range_x):
    """
    Calcula os possíveis ângulos de lançamento (theta).
    Retorna: (angle1, angle2, steps_string)
    """
    steps = f"--- Passo-a-passo: Ângulo de Lançamento ---\n"
    steps += f"1. Dados Fornecidos:\n"
    steps += f"   Velocidade Inicial (v0): {v0} m/s\n"
    steps += f"   Altura Inicial (y0): {y0} m\n"
    steps += f"   Altura Final (y_final): {y_final} m\n"
    steps += f"   Alcance (x): {range_x} m\n"

    if v0 == 0:
        raise ValueError("Velocidade inicial (v0) não pode ser zero.")

    steps += f"\n2. Usando a Eq. da Trajetória (com 1/cos²θ = 1 + tan²θ):\n"
    steps += f"   y = y0 + x*tan(θ) - (g*x²) / (2*v0²)*(1 + tan²(θ))\n"
    steps += f"   Para simplificar, definimos u = tan(θ) e K = (g*x²) / (2*v0²)\n"

    K = (GRAVITY * range_x**2) / (2 * v0**2)
    steps += f"   K = ({GRAVITY} * {range_x}²) / (2 * {v0}²) = {K:.2f}\n"
    steps += f"   Substituindo u e K:\n"
    steps += f"   {y_final} = {y0} + {range_x}*u - {K:.2f}*(1 + u²)\n"
    steps += f"   Reorganizando para a forma a*u² + b*u + c = 0:\n"

    a = K
    b = -range_x
    c = y_final - y0 + K
    steps += f"   ({a:.2f})u² + ({b:.2f})u + ({c:.2f}) = 0\n"

    steps += f"3. Resolvendo Bhaskara para u (u = tan(θ)):\n"
    discriminant = b**2 - 4*a*c
    steps += f"   Δ = ({b:.2f})² - 4*({a:.2f})*({c:.2f}) = {discriminant:.2f}\n"

    if discriminant < 0:
        raise ValueError("Alvo fora de alcance para esta velocidade inicial.")

    u1 = (-b + math.sqrt(discriminant)) / (2*a)
    u2 = (-b - math.sqrt(discriminant)) / (2*a)
    steps += f"   u1 = ({-b:.2f} + √{discriminant:.2f}) / (2*{a:.2f}) = {u1:.2f}\n"
    steps += f"   u2 = ({-b:.2f} - √{discriminant:.2f}) / (2*{a:.2f}) = {u2:.2f}\n"

    steps += f"\n4. Encontrando os Ângulos (θ = atan(u)):\n"
    angle1 = math.degrees(math.atan(u1))
    angle2 = math.degrees(math.atan(u2))
    steps += f"   θ1 = atan({u1:.2f}) = {angle1:.2f}°\n"
    steps += f"   θ2 = atan({u2:.2f}) = {angle2:.2f}°\n"

    steps += f"\n--- RESULTADO FINAL ---\n"
    steps += f"   Ângulos Possíveis: {angle1:.2f}° ou {angle2:.2f}°\n"

    return angle1, angle2, steps


def calculate_from_vertex(y0, Sx, Sy):
    """
    Calcula v0 e theta a partir do vértice (Sx, Sy) e altura inicial y0.
    Retorna: (v0, theta_deg, steps_string)
    """
    steps = f"--- Passo-a-passo: Cálculo pelo Vértice ---\n"
    steps += f"1. Dados Fornecidos:\n"
    steps += f"   Altura Inicial (y0): {y0} m\n"
    steps += f"   Posição X do Vértice (Sx): {Sx} m\n"
    steps += f"   Posição Y do Vértice (Sy, Altura Máx.): {Sy} m\n"

    if Sy < y0:
        raise ValueError(
            "Altura do vértice (Sy) deve ser >= altura inicial (y0).")
    if Sx < 0:
        raise ValueError("Alcance até o vértice (Sx) deve ser >= 0.")

    steps += f"\n2. Encontrando v0y (Velocidade Vertical Inicial):\n"
    steps += f"   No vértice, vy = 0. Usamos Torricelli:\n"
    steps += f"   vy² = v0y² - 2*g*(Sy - y0)\n"
    steps += f"   0 = v0y² - 2*({GRAVITY})*({Sy} - {y0})\n"
    v0y = math.sqrt(2 * GRAVITY * (Sy - y0))
    steps += f"   v0y = √[2 * {GRAVITY} * ({Sy - y0:.2f})] = {v0y:.2f} m/s\n"

    steps += f"\n3. Encontrando o Tempo até o Pico (t_pico):\n"
    steps += f"   Usamos a Equação da Velocidade (vy = 0):\n"
    steps += f"   vy = v0y - g*t\n"
    steps += f"   0 = {v0y:.2f} - {GRAVITY} * t_pico\n"

    if v0y == 0:  # Caso de lançamento horizontal (Sy = y0)
        if Sx != 0:
            raise ValueError(
                "Se Sy=y0 (vértice na altura inicial), Sx deve ser 0.")
        t_peak = 0.0
    else:
        t_peak = v0y / GRAVITY

    steps += f"   t_pico = {v0y:.2f} / {GRAVITY} = {t_peak:.2f} s\n"

    steps += f"\n4. Encontrando v0x (Velocidade Horizontal):\n"
    steps += f"   O movimento horizontal é constante: Sx = v0x * t_pico\n"
    if t_peak == 0:
        if Sx == 0:
            # Lançamento horizontal (v0y=0) mas sem movimento (Sx=0)
            # Implica v0x = 0 também.
            v0x = 0.0
            steps += f"   Como t_pico é 0, v0x também é 0.\n"
        else:
            raise ValueError("Contradição: Tempo de pico é 0 mas Sx > 0.")
    else:
        v0x = Sx / t_peak
        steps += f"   v0x = Sx / t_pico = {Sx} / {t_peak:.2f} = {v0x:.2f} m/s\n"

    steps += f"\n5. Combinando v0x e v0y para encontrar v0 e θ:\n"
    v0 = math.sqrt(v0x**2 + v0y**2)
    steps += f"   v0 = √(v0x² + v0y²) = √({v0x:.2f}² + {v0y:.2f}²) = {v0:.2f} m/s\n"

    # atan2 é mais robusto para todos os quadrantes, evita divisão por zero
    theta_rad = math.atan2(v0y, v0x)
    theta_deg = math.degrees(theta_rad)
    steps += f"   θ = atan(v0y / v0x) = atan({v0y:.2f} / {v0x:.2f}) = {theta_deg:.2f}°\n"

    steps += f"\n--- RESULTADO FINAL (Dados de Lançamento) ---\n"
    steps += f"   Velocidade Inicial (v0): {v0:.2f} m/s\n"
    steps += f"   Ângulo de Lançamento (θ): {theta_deg:.2f}°\n"

    return v0, theta_deg, steps


# --- Classe da Aplicação GUI ---

class PhysicsCalculatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Calculadora de Lançamento (Física 1)")
        self.root.geometry("900x700")  # Janela um pouco mais larga

        self.calculation_mode = tk.StringVar()
        self.entries = {}
        self.labels = {}

        # --- Criação dos Widgets ---
        self.create_mode_selector()

        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill='both', expand=True, side='left', anchor='n')

        self.create_input_frame(main_frame)
        self.create_action_frame(main_frame)
        self.create_result_frame(main_frame)  # <-- MODIFICADO

        self.plot_frame = ttk.LabelFrame(
            self.root, text="Gráfico da Trajetória", padding="10")
        self.plot_frame.pack(fill='both', expand=True,
                             side='right', padx=10, pady=5)

        self.fig, self.ax = plt.subplots(figsize=(5, 5))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill='both', expand=True)
        self.clear_plot()

        self.calculation_mode.set("Calcular Alcance e Tempo de Voo")
        self.update_gui_state(None)

    def create_mode_selector(self):
        """Cria o menu para selecionar o modo de cálculo."""
        frame = ttk.Frame(self.root, padding="10")
        frame.pack(fill='x')

        label = ttk.Label(frame, text="O que você quer calcular?")
        label.pack(side='left', padx=5)

        modes = [
            "Calcular Alcance e Tempo de Voo",
            "Calcular Altura Máxima",
            "Calcular Velocidade Inicial",
            "Calcular Ângulo de Lançamento",
            "Calcular v0 e Ângulo (pelo Vértice)"
        ]

        mode_menu = ttk.Combobox(
            frame, textvariable=self.calculation_mode, values=modes, state='readonly')
        mode_menu.pack(side='left', fill='x', expand=True)
        mode_menu.bind('<<ComboboxSelected>>', self.update_gui_state)

    def create_input_frame(self, parent_frame):
        """Cria os campos de entrada para todas as variáveis."""
        frame = ttk.LabelFrame(parent_frame, text="Informações", padding="10")
        frame.pack(fill='x', padx=10, pady=5)

        input_vars = {
            "v0": "Velocidade Inicial (v0) [m/s]:",
            "theta": "Ângulo (θ) [graus]:",
            "y0": "Altura Inicial (y0) [m]:",
            "y_final": "Altura Final (y) [m]:",
            "range_x": "Alcance (x) [m]:",
            "Sx": "Alcance até o Vértice (Sx) [m]:",
            "Sy": "Altura do Vértice (Sy) [m]:"
        }

        for i, (key, text) in enumerate(input_vars.items()):
            label = ttk.Label(frame, text=text)
            label.grid(row=i, column=0, sticky='w', padx=5, pady=5)

            entry = ttk.Entry(frame)
            entry.grid(row=i, column=1, sticky='ew', padx=5, pady=5)

            self.labels[key] = label
            self.entries[key] = entry

        frame.columnconfigure(1, weight=1)

    def create_action_frame(self, parent_frame):
        """Cria o botão de calcular."""
        frame = ttk.Frame(parent_frame, padding="10")
        frame.pack(fill='x')

        calculate_button = ttk.Button(
            frame, text="Calcular", command=self.perform_calculation)
        calculate_button.pack(expand=True)

    def create_result_frame(self, parent_frame):
        """Cria a área para mostrar os resultados (agora com Text e Scrollbar)."""
        frame = ttk.LabelFrame(
            parent_frame, text="Passo a Passo e Resultados", padding="10")
        frame.pack(fill='both', expand=True, padx=10, pady=10)

        # Cria a barra de rolagem
        scrollbar = ttk.Scrollbar(frame, orient='vertical')

        # Cria o widget de Texto
        self.result_text_widget = tk.Text(frame,
                                          wrap='word',
                                          height=15,
                                          # Fonte monoespaçada
                                          font=("Courier New", 10),
                                          yscrollcommand=scrollbar.set)

        # Configura a barra de rolagem para controlar o texto
        scrollbar.config(command=self.result_text_widget.yview)

        # Empacota os widgets
        scrollbar.pack(side='right', fill='y')
        self.result_text_widget.pack(fill='both', expand=True)

        # Inicia com um texto e desabilitado para edição
        self.set_result_text("Preencha os dados e clique em 'Calcular'.")

    def set_result_text(self, text):
        """Limpa e insere texto no widget de resultado."""
        self.result_text_widget.config(state='normal')  # Habilita para escrita
        self.result_text_widget.delete('1.0', 'end')
        self.result_text_widget.insert('1.0', text)
        # Desabilita para o usuário não editar
        self.result_text_widget.config(state='disabled')

    def update_gui_state(self, event):
        """Habilita/Desabilita campos de entrada com base no modo selecionado."""
        mode = self.calculation_mode.get()

        required_inputs = {
            "Calcular Alcance e Tempo de Voo": ["v0", "theta", "y0", "y_final"],
            "Calcular Altura Máxima": ["v0", "theta", "y0"],
            "Calcular Velocidade Inicial": ["theta", "y0", "y_final", "range_x"],
            "Calcular Ângulo de Lançamento": ["v0", "y0", "y_final", "range_x"],
            "Calcular v0 e Ângulo (pelo Vértice)": ["y0", "y_final", "Sx", "Sy"]
        }

        needed = required_inputs.get(mode, [])

        for key, entry in self.entries.items():
            label = self.labels[key]
            if key in needed:
                entry.config(state='normal')
                label.config(state='normal')
            else:
                entry.delete(0, 'end')
                entry.config(state='disabled')
                label.config(state='disabled')

        self.set_result_text("Preencha os dados e clique em 'Calcular'.")
        self.clear_plot()

    def get_float_from_entry(self, key, default_value=None):
        """Função auxiliar para ler e converter o valor de uma entry."""
        label_text = self.labels[key].cget("text")

        try:
            value_str = self.entries[key].get()
            if not value_str:
                if default_value is not None:
                    self.entries[key].insert(0, str(default_value))
                    return default_value
                raise ValueError(f"O campo '{label_text}' está vazio.")
            # Aceita vírgula ou ponto
            return float(value_str.replace(',', '.'))
        except ValueError:
            raise ValueError(
                f"Entrada inválida para '{label_text}'. Use apenas números.")

    def plot_trajectory(self, v0, theta_deg, y0, time_of_flight):
        """Plota a trajetória do projétil."""
        self.ax.clear()

        theta_rad = math.radians(theta_deg)
        v0x = v0 * math.cos(theta_rad)
        v0y = v0 * math.sin(theta_rad)

        if not isinstance(time_of_flight, (int, float)) or time_of_flight <= 0:
            self.ax.text(0.5, 0.5, 'Gráfico indisponível (tempo de voo inválido)',
                         horizontalalignment='center', verticalalignment='center',
                         transform=self.ax.transAxes)
            self.canvas.draw()
            return

        t = np.linspace(0, time_of_flight, 100)
        x = v0x * t
        y = y0 + v0y * t - 0.5 * GRAVITY * t**2

        self.ax.plot(x, y, label='Trajetória')
        self.ax.plot(0, y0, 'go', label=f'Lançamento (0, {y0})')
        self.ax.plot(x[-1], y[-1], 'ro',
                     label=f'Aterrissagem ({x[-1]:.1f}, {y[-1]:.1f})')

        if v0y > 0:
            t_peak = v0y / GRAVITY
            if t_peak <= time_of_flight:
                x_peak = v0x * t_peak
                y_peak = y0 + v0y * t_peak - 0.5 * GRAVITY * t_peak**2
                self.ax.plot(x_peak, y_peak, 'b*', markersize=8,
                             label=f'Vértice ({x_peak:.1f}, {y_peak:.1f})')

        self.ax.set_xlabel("Distância Horizontal (m)")
        self.ax.set_ylabel("Altura (m)")
        self.ax.set_title("Trajetória do Projétil")
        self.ax.grid(True)
        self.ax.axhline(0, color='black', linewidth=0.5)

        # Ajusta limites para visualização
        min_y = min(np.min(y), 0, y0)
        max_y = max(np.max(y), y0)
        self.ax.set_ylim(bottom=min_y - (max_y - min_y) * 0.1,  # 10% de margem
                         top=max_y + (max_y - min_y) * 0.1)

        min_x = min(np.min(x), 0)
        max_x = max(np.max(x), 0.1)  # Evita x=0
        self.ax.set_xlim(left=min_x - (max_x - min_x) * 0.1,
                         right=max_x + (max_x - min_x) * 0.1)

        self.ax.legend()
        self.canvas.draw()

    def clear_plot(self):
        """Limpa o gráfico."""
        self.ax.clear()
        self.ax.set_xlabel("Distância Horizontal (m)")
        self.ax.set_ylabel("Altura (m)")
        self.ax.set_title("Gráfico da Trajetória")
        self.ax.grid(True)
        self.ax.axhline(0, color='black', linewidth=0.5)
        self.canvas.draw()

    def perform_calculation(self):
        """Orquestra o processo de cálculo e exibe o resultado."""
        mode = self.calculation_mode.get()
        self.clear_plot()

        try:
            if mode == "Calcular Alcance e Tempo de Voo":
                v0 = self.get_float_from_entry("v0")
                theta = self.get_float_from_entry("theta")
                y0 = self.get_float_from_entry("y0")
                y_final = self.get_float_from_entry(
                    "y_final", default_value=0.0)

                range_x, time, steps = calculate_range_and_time(
                    v0, theta, y0, y_final)
                self.set_result_text(steps)
                self.plot_trajectory(v0, theta, y0, time)

            elif mode == "Calcular Altura Máxima":
                v0 = self.get_float_from_entry("v0")
                theta = self.get_float_from_entry("theta")
                y0 = self.get_float_from_entry("y0")

                max_h, time_peak, steps = calculate_max_height(v0, theta, y0)
                self.set_result_text(steps)

                try:
                    # Tenta calcular o tempo total para plotar o gráfico completo
                    _, time_total, _ = calculate_range_and_time(
                        v0, theta, y0, 0.0)  # Assume queda até y=0
                    self.plot_trajectory(v0, theta, y0, time_total)
                except Exception:
                    # Se falhar (ex: lançamento para baixo), plota só até o pico
                    self.plot_trajectory(v0, theta, y0, time_peak)

            elif mode == "Calcular Velocidade Inicial":
                theta = self.get_float_from_entry("theta")
                y0 = self.get_float_from_entry("y0")
                y_final = self.get_float_from_entry(
                    "y_final", default_value=0.0)
                range_x = self.get_float_from_entry("range_x")

                v0, steps = calculate_initial_velocity(
                    theta, y0, y_final, range_x)
                self.set_result_text(steps)

                _, time_flight, _ = calculate_range_and_time(
                    v0, theta, y0, y_final)
                self.plot_trajectory(v0, theta, y0, time_flight)

            elif mode == "Calcular Ângulo de Lançamento":
                v0 = self.get_float_from_entry("v0")
                y0 = self.get_float_from_entry("y0")
                y_final = self.get_float_from_entry(
                    "y_final", default_value=0.0)
                range_x = self.get_float_from_entry("range_x")

                angle1, angle2, steps = calculate_angle(
                    v0, y0, y_final, range_x)
                self.set_result_text(steps)
                self.plot_two_trajectories(
                    v0, angle1, y0, y_final, v0, angle2, y0, y_final, range_x)

            elif mode == "Calcular v0 e Ângulo (pelo Vértice)":
                y0 = self.get_float_from_entry("y0")
                y_final = self.get_float_from_entry(
                    "y_final", default_value=0.0)
                Sx = self.get_float_from_entry("Sx")
                Sy = self.get_float_from_entry("Sy")

                v0, theta, steps_vertex = calculate_from_vertex(y0, Sx, Sy)

                # Agora, calcula o alcance total e tempo total para o y_final fornecido
                range_x, time, steps_total = calculate_range_and_time(
                    v0, theta, y0, y_final)

                # Remove o cabeçalho do segundo passo-a-passo para evitar redundância
                steps_total_cleaned = "\n".join(steps_total.split("\n")[1:])

                final_explanation = steps_vertex
                final_explanation += f"\n\n--- Com v0={v0:.2f} m/s e θ={theta:.2f}° encontrados, calculamos a trajetória completa até y_final={y_final} m ---\n"
                final_explanation += steps_total_cleaned

                self.set_result_text(final_explanation)
                self.plot_trajectory(v0, theta, y0, time)

        except ValueError as ve:
            messagebox.showerror("Erro de Entrada ou Cálculo", str(ve))
            self.set_result_text(f"Erro: {ve}")
        except Exception as e:
            messagebox.showerror("Erro Inesperado", f"Ocorreu um erro: {e}")
            self.set_result_text(f"Erro: {e}")

    def plot_two_trajectories(self, v0_a, theta_deg_a, y0_a, y_final_a, v0_b, theta_deg_b, y0_b, y_final_b, target_range_x):
        """Plota duas trajetórias, usado para o cálculo de ângulo."""
        self.ax.clear()

        all_y = [y0_a]
        all_x = [0, target_range_x]

        # Plot Trajetória 1
        try:
            _, time_a, _ = calculate_range_and_time(
                v0_a, theta_deg_a, y0_a, y_final_a)
            theta_rad_a = math.radians(theta_deg_a)
            v0x_a = v0_a * math.cos(theta_rad_a)
            v0y_a = v0_a * math.sin(theta_rad_a)
            t_a = np.linspace(0, time_a, 100)
            x_a = v0x_a * t_a
            y_a = y0_a + v0y_a * t_a - 0.5 * GRAVITY * t_a**2
            self.ax.plot(
                x_a, y_a, label=f'Trajetória 1 (θ={theta_deg_a:.2f}°)')
            all_y.extend(y_a)
            all_x.extend(x_a)
        except Exception:
            pass  # Ignora se um dos ângulos for impossível

        # Plot Trajetória 2
        try:
            _, time_b, _ = calculate_range_and_time(
                v0_b, theta_deg_b, y0_b, y_final_b)
            theta_rad_b = math.radians(theta_deg_b)
            v0x_b = v0_b * math.cos(theta_rad_b)
            v0y_b = v0_b * math.sin(theta_rad_b)
            t_b = np.linspace(0, time_b, 100)
            x_b = v0x_b * t_b
            y_b = y0_b + v0y_b * t_b - 0.5 * GRAVITY * t_b**2
            self.ax.plot(
                x_b, y_b, label=f'Trajetória 2 (θ={theta_deg_b:.2f}°)', linestyle='--')
            all_y.extend(y_b)
            all_x.extend(x_b)
        except Exception:
            pass

        self.ax.plot(0, y0_a, 'go', label='Ponto de Lançamento')
        self.ax.plot(target_range_x, y_final_a, 'ro',
                     label='Ponto Final Desejado')

        self.ax.set_xlabel("Distância Horizontal (m)")
        self.ax.set_ylabel("Altura (m)")
        self.ax.set_title("Trajetórias para Ângulos Possíveis")
        self.ax.grid(True)
        self.ax.axhline(0, color='black', linewidth=0.5)

        # Ajusta limites
        min_y = min(np.min(all_y), 0)
        max_y = max(np.max(all_y), 0.1)
        self.ax.set_ylim(bottom=min_y - (max_y - min_y) * 0.1,
                         top=max_y + (max_y - min_y) * 0.1)

        min_x = min(np.min(all_x), 0)
        max_x = max(np.max(all_x), 0.1)
        self.ax.set_xlim(left=min_x - (max_x - min_x) * 0.1,
                         right=max_x + (max_x - min_x) * 0.1)

        self.ax.legend()
        self.canvas.draw()


# --- Ponto de Entrada Principal ---
if __name__ == "__main__":
    main_window = tk.Tk()
    app = PhysicsCalculatorApp(main_window)
    main_window.mainloop()
