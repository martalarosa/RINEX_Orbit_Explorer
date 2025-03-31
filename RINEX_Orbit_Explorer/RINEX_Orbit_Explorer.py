# данный код рассчитывает и строит эволюции орбиты и параметров спутника, используя данные из навигационного послания в RINEX-формате

import math
import pandas as pd
import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# константы
EARTH_GRAVITY_PARAM = 3.986005 * 10**14  # гравитационный параметр Земли (м^3/с^2)
PI = math.pi
EARTH_ROTATION_RATE = 7.2921151467e-5  # скорость вращения Земли (рад/с)

# функция для чтения RINEX-файла и извлечения данных для спутника
def parse_rinex_file(file_path, satellite_id):
    with open(file_path, 'r', encoding='latin-1') as file:
        lines = file.readlines()
    
    satellite_records = []  # список для хранения данных спутника
    is_collecting = False  # флаг для начала сбора данных
    for line in lines:
        if line.startswith(satellite_id): # если строка начинается с ID спутника, начинаем сбор данных
            is_collecting = True
        elif is_collecting and line.startswith('G'): # если встречаем следующий спутник, прекращаем сбор
            break
        if is_collecting:
            formatted_line = re.sub(r'(\d)(-)', r'\1 -', line.strip())  # исправляем формат чисел
            satellite_records.append(formatted_line)
    
    return satellite_records

# функция для извлечения параметров орбиты
def extract_orbital_params(satellite_records):
    try:
        mean_motion_delta = float(satellite_records[1].split()[2].replace('D', 'E'))  # изменение среднего движения спутника
        mean_anomaly_at_epoch = float(satellite_records[1].split()[3].replace('D', 'E'))  # средняя аномалия в момент эпохи
        eccentricity = float(satellite_records[2].split()[1].replace('D', 'E'))  # эксцентриситет орбиты
        semi_major_axis_sqrt = float(satellite_records[2].split()[3].replace('D', 'E'))  # квадратный корень из большой полуоси
        time_of_ephemeris = float(satellite_records[3].split()[0].replace('D', 'E'))  # время эпохи эфемерид
        perigee_argument = float(satellite_records[4].split()[2].replace('D', 'E'))  # аргумент перигея
        right_ascension_at_epoch = float(satellite_records[3].split()[2].replace('D', 'E'))  # прямое восхождение (долгота) восходящего узла в момент эпохи
        right_ascension_rate = float(satellite_records[4].split()[3].replace('D', 'E'))  # скорость изменения прямого восхождения
        inclination_at_epoch = float(satellite_records[4].split()[0].replace('D', 'E'))  # наклонение орбиты в момент эпохи
        inclination_rate = float(satellite_records[5].split()[0].replace('D', 'E'))  # скорость изменения наклонения орбиты

        # корректирующие поправки к параметрам орбиты
        correction_params = {
            'C_uc': float(satellite_records[2].split()[0].replace('D', 'E')),  # косинусная поправка к аргументу широты
            'C_us': float(satellite_records[2].split()[2].replace('D', 'E')),  # синусная поправка к аргументу широты
            'C_rc': float(satellite_records[4].split()[1].replace('D', 'E')),  # косинусная поправка к радиусу орбиты
            'C_rs': float(satellite_records[1].split()[1].replace('D', 'E')),  # синусная поправка к радиусу орбиты
            'C_ic': float(satellite_records[3].split()[1].replace('D', 'E')),  # косинусная поправка к наклонению орбиты
            'C_is': float(satellite_records[3].split()[3].replace('D', 'E'))   # синусная поправка к наклонению орбиты
        }
    except IndexError:
        raise ValueError("Ошибка: Неверный формат данных в файле.")
    
    params = (mean_motion_delta, mean_anomaly_at_epoch, eccentricity, semi_major_axis_sqrt, 
              time_of_ephemeris, perigee_argument, right_ascension_at_epoch, right_ascension_rate, 
              inclination_at_epoch, inclination_rate, correction_params)
    
    print("Извлеченные параметры орбиты:")
    param_names = [
        "Изменение среднего движения спутника", "Средняя аномалия в момент эпохи", "Эксцентриситет орбиты", "Квадратный корень из большой полуоси", 
        "Время эпохи эфемерид", "Аргумент перигея", "Прямое восхождение восходящего узла в момент эпохи", "Скорость изменения прямого восхождения", 
        "Наклонение орбиты в момент эпохи", "Скорость изменения наклонения орбиты", "Корректирующие поправки"
    ]
    for param, name in zip(params, param_names):
        print(f"{name}: {param}")
    
    return params

# функция для вычисления орбиты
def calculate_orbit(mean_motion_delta, mean_anomaly_at_epoch, eccentricity, semi_major_axis_sqrt, 
                     time_of_ephemeris, perigee_argument, right_ascension_at_epoch, right_ascension_rate, 
                     inclination_at_epoch, inclination_rate, correction_params):
    semi_major_axis = semi_major_axis_sqrt ** 2 # вычисляем большую полуось
    mean_motion = math.sqrt(EARTH_GRAVITY_PARAM / semi_major_axis**3) + mean_motion_delta # среднее движение (скорость движения спутника по орбите)
    orbit_duration = 24 * 60 * 60  # период обращения спутника
    time_step = 60  # шаг по времени в секундах

    x_positions, y_positions, z_positions = [], [], []
    orbital_params = {
        't': [], 'M': [], 'E': [], 'true_anomaly': [], 'u': [], 'r': [],
        'i': [], 'omega': [], 'delta_u': [], 'delta_r': [], 'delta_i': []
    }

   # запускаем цикл, чтобы вычислить траекторию
    for time in range(int(time_of_ephemeris), int(time_of_ephemeris + orbit_duration), time_step):
        # решаем уравнение Кеплера для эксцентрической аномалии
        M = mean_anomaly_at_epoch + mean_motion * (time - time_of_ephemeris)
        E = M
        for _ in range(100):
            E = M + eccentricity * math.sin(E)

        # находим истинную аномалию (положение спутника в орбите)
        true_anomaly = 2 * math.atan(math.sqrt((1 + eccentricity) / (1 - eccentricity)) * math.tan(E / 2))
        if true_anomaly < 0:
            true_anomaly += 2 * math.pi
        
        # корректируем орбитальные параметры с учётом возмущений
        argument_of_latitude = true_anomaly + perigee_argument
        delta_u = correction_params['C_us'] * math.sin(2 * argument_of_latitude) + correction_params['C_uc'] * math.cos(2 * argument_of_latitude)
        delta_r = correction_params['C_rs'] * math.sin(2 * argument_of_latitude) + correction_params['C_rc'] * math.cos(2 * argument_of_latitude)
        delta_i = correction_params['C_is'] * math.sin(2 * argument_of_latitude) + correction_params['C_ic'] * math.cos(2 * argument_of_latitude)

        u = argument_of_latitude + delta_u
        r = semi_major_axis * (1 - eccentricity * math.cos(E)) + delta_r
        i = inclination_at_epoch + delta_i + inclination_rate * (time - time_of_ephemeris)
        omega = right_ascension_at_epoch + (right_ascension_rate - EARTH_ROTATION_RATE) * (time - time_of_ephemeris) - EARTH_ROTATION_RATE * time_of_ephemeris

        # переводим в прямоугольные координаты ECEF (геоцентрическая СК)
        x_orbit = r * math.cos(u)
        y_orbit = r * math.sin(u)

        X = x_orbit * math.cos(omega) - y_orbit * math.cos(i) * math.sin(omega)
        Y = x_orbit * math.sin(omega) + y_orbit * math.cos(i) * math.cos(omega)
        Z = y_orbit * math.sin(i)

        x_positions.append(X)
        y_positions.append(Y)
        z_positions.append(Z)

        # сохраняем параметры
        orbital_params['t'].append(time - time_of_ephemeris)
        orbital_params['M'].append(M)
        orbital_params['E'].append(E)
        orbital_params['true_anomaly'].append(true_anomaly)
        orbital_params['u'].append(u)
        orbital_params['r'].append(r)
        orbital_params['i'].append(i)
        orbital_params['omega'].append(omega)
        orbital_params['delta_u'].append(delta_u)
        orbital_params['delta_r'].append(delta_r)
        orbital_params['delta_i'].append(delta_i)

    # получаем списки x_positions, y_positions, z_positions, содержащие траекторию спутника
    return x_positions, y_positions, z_positions, true_anomaly, orbital_params

# функция для вывода истинной аномалии
def display_true_anomaly(true_anomaly):
    print(f"Истинная аномалия в момент эпохи:", true_anomaly)

# визуализация эволюции орбиты (ECEF) в 3D- и 2D-плоскостях
def plot_satellite_orbit(x_positions, y_positions, z_positions):
    fig = plt.figure(figsize=(15, 5))
    
    ax1 = fig.add_subplot(221, projection='3d')
    ax1.plot(x_positions, y_positions, z_positions, color='purple', label='Орбита КА')
    ax1.scatter(0, 0, 0, color='blue', label='Земля', s=88)
    ax1.set_title('Эволюция координат (орбиты) в трёхмерной плоскости')
    ax1.legend()
    
    ax2 = fig.add_subplot(222)
    ax2.plot(x_positions, y_positions, color='red')
    ax2.set_title('Эволюция координат (орбиты) в плоскости X-Y')
    
    ax3 = fig.add_subplot(223)
    ax3.plot(x_positions, z_positions, color='yellow')
    ax3.set_title('Эволюция координат (орбиты) в плоскости X-Z')
    
    ax4 = fig.add_subplot(224)
    ax4.plot(y_positions, z_positions, color='orange')
    ax4.set_title('Эволюция координат (орбиты) в плоскости Y-Z')
    
    plt.show()

# визуализация эволюций орбитальных параметров
def plot_orbital_parameters(params):
    t = [ti / 3600 for ti in params['t']]

    fig, axs = plt.subplots(5, 2, figsize=(14, 12))
    fig.suptitle("Эволюции орбитальных параметров", fontsize=16)

    axs[0, 0].plot(t, [math.degrees(val) for val in params['M']], label='M')
    axs[0, 0].set_title("Средняя аномалия (°)")

    axs[0, 1].plot(t, [math.degrees(val) for val in params['E']], label='E')
    axs[0, 1].set_title("Эксцентрическая аномалия (°)")

    axs[1, 0].plot(t, [math.degrees(val) for val in params['true_anomaly']], label='ν')
    axs[1, 0].set_title("Истинная аномалия (°)")

    axs[1, 1].plot(t, [math.degrees(val) for val in params['u']], label='u')
    axs[1, 1].set_title("Аргумент широты (°)")

    axs[2, 0].plot(t, params['r'], label='r')
    axs[2, 0].set_title("Радиус орбиты (м)")

    axs[2, 1].plot(t, [math.degrees(val) for val in params['i']], label='i')
    axs[2, 1].set_title("Наклонение (°)")

    axs[3, 0].plot(t, [math.degrees(val) for val in params['omega']], label='Ω')
    axs[3, 0].set_title("Долгота восходящего узла (°)")

    axs[3, 1].plot(t, [math.degrees(val) for val in params['delta_u']], label='δu')
    axs[3, 1].set_title("Поправка δu (°)")

    axs[4, 0].plot(t, params['delta_r'], label='δr')
    axs[4, 0].set_title("Поправка δr (м)")

    axs[4, 1].plot(t, [math.degrees(val) for val in params['delta_i']], label='δi')
    axs[4, 1].set_title("Поправка δi (°)")

    for ax_row in axs:
        for ax in ax_row:
            ax.set_xlabel("Время (ч)")
            ax.grid(True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# функция для вывода координат (числовых данных траектории) в таблицу в консоль и записи в файл 
def display_orbit_coordinates(x_positions, y_positions, z_positions, params, filename="orbit_coordinates.csv"):
    # определяем начальные координаты
    X_epoch = x_positions[0]
    Y_epoch = y_positions[0]
    Z_epoch = z_positions[0]
    
    # создаём словарь с данными
    data = {
        "Время (с)": [params[4]] + list(range(int(params[4]), int(params[4]) + len(x_positions) * 60, 60)),
        "Координата X (м)": [X_epoch] + x_positions,
        "Координата Y (м)": [Y_epoch] + y_positions,
        "Координата Z (м)": [Z_epoch] + z_positions
    }

    df = pd.DataFrame(data) # преобразуем данные в DataFrame

    df.to_csv(filename, index=False)  # запишем данные в CSV файл без индексов

    print("\nПолный список координат спутника в системе ECEF:")
    print(df.to_string(index=False))  # печать в консоль без индексов

    print(f"\nДанные записаны в файл {filename}.")


# основной код для обработки спутниковых данных
if __name__ == '__main__':
    rinex_file_path = 'C:/Users/МАРТА/RINEX_Orbit_Explorer/RINEX_Orbit_Explorer/Brdc0020.25n'
    satellite_id = 'G01' # ID спутника, данные которого необходимо обработать

    satellite_records = parse_rinex_file(rinex_file_path, satellite_id) # чтение и извлечение данных
    params = extract_orbital_params(satellite_records) # извлечение параметров орбиты

    x, y, z, true_anomaly, orbital_params = calculate_orbit(*params) # вычисление орбиты и параметров
    display_true_anomaly(true_anomaly) # вывод истинной аномалии

    display_orbit_coordinates(x, y, z, params, filename="orbit_coordinates.csv") # вывод координат
    plot_satellite_orbit(x, y, z) # визуализация эволюции орбиты спутника
    plot_orbital_parameters(orbital_params) # визуализация эволюций орбитальных параметров