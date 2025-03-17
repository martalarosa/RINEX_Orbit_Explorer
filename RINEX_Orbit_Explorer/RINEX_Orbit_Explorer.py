# данный код рассчитывает и строит траекторию движения спутника, используя данные из навигационного послания в RINEX-формате

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
        right_ascension_at_epoch = float(satellite_records[3].split()[2].replace('D', 'E'))  # прямое восхождение восходящего узла в момент эпохи
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
    
    print("\nИзвлеченные параметры орбиты:")
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
    semi_major_axis = semi_major_axis_sqrt ** 2 # вычисляем большую полуось орбиты
    mean_motion = math.sqrt(EARTH_GRAVITY_PARAM / semi_major_axis**3) + mean_motion_delta # среднее движение (скорость движения спутника по орбите)
    orbit_duration = 12 * 60 * 60  # период обращения спутника (12 часов)
    time_step = 60  # шаг по времени в секундах

    x_positions, y_positions, z_positions = [], [], []
    
   # запускаем цикл, чтобы вычислить траекторию на 12 часов
    for time in range(int(time_of_ephemeris), int(time_of_ephemeris + orbit_duration), time_step):
        # решаем уравнение Кеплера для эксцентрической аномалии
        mean_anomaly = mean_anomaly_at_epoch + mean_motion * (time - time_of_ephemeris)
        eccentric_anomaly = mean_anomaly
        for _ in range(3):
            eccentric_anomaly = mean_anomaly + eccentricity * math.sin(eccentric_anomaly)
        # находим истинную аномалию (положение спутника в орбите)
        true_anomaly = 2 * math.atan(math.sqrt((1 + eccentricity) / (1 - eccentricity)) * math.tan(eccentric_anomaly / 2))
        # корректируем орбитальные параметры с учётом возмущений
        argument_of_latitude = true_anomaly + perigee_argument
        delta_u = correction_params['C_us'] * math.sin(2 * argument_of_latitude) + correction_params['C_uc'] * math.cos(2 * argument_of_latitude)
        delta_r = correction_params['C_rs'] * math.sin(2 * argument_of_latitude) + correction_params['C_rc'] * math.cos(2 * argument_of_latitude)
        delta_i = correction_params['C_is'] * math.sin(2 * argument_of_latitude) + correction_params['C_ic'] * math.cos(2 * argument_of_latitude)
        
        corrected_latitude = argument_of_latitude + delta_u
        orbit_radius = semi_major_axis * (1 - eccentricity * math.cos(eccentric_anomaly)) + delta_r
        inclination = inclination_at_epoch + delta_i + inclination_rate * (time - time_of_ephemeris)
        # переводим в прямоугольные координаты ECEF
        x_orbit = orbit_radius * math.cos(corrected_latitude)
        y_orbit = orbit_radius * math.sin(corrected_latitude)
        
        right_ascension = right_ascension_at_epoch + (right_ascension_rate - EARTH_ROTATION_RATE) * (time - time_of_ephemeris) - EARTH_ROTATION_RATE * time_of_ephemeris
        X = x_orbit * math.cos(right_ascension) - y_orbit * math.cos(inclination) * math.sin(right_ascension)
        Y = x_orbit * math.sin(right_ascension) + y_orbit * math.cos(inclination) * math.cos(right_ascension)
        Z = y_orbit * math.sin(inclination)
        
        x_positions.append(X)
        y_positions.append(Y)
        z_positions.append(Z)
    # получаем списки x_positions, y_positions, z_positions, содержащие траекторию спутника
    return x_positions, y_positions, z_positions

# визуализация орбиты в 3D- и 2D-плоскостях
def plot_satellite_orbit(x_positions, y_positions, z_positions):
    fig = plt.figure(figsize=(15, 5))
    
    ax1 = fig.add_subplot(131, projection='3d')
    ax1.plot(x_positions, y_positions, z_positions, color='red', label='Орбита КА')
    ax1.scatter(0, 0, 0, color='yellow', label='Земля', s=100)
    ax1.set_title('Трёхмерная плоскость')
    ax1.legend()
    
    ax2 = fig.add_subplot(132)
    ax2.plot(x_positions, y_positions, color='red')
    ax2.set_title('XY плоскость')
    
    ax3 = fig.add_subplot(133)
    ax3.plot(x_positions, z_positions, color='red')
    ax3.set_title('XZ плоскость')
    
    plt.show()

# функция для вывода координат в таблице (числовые данные траектории)
def display_orbit_coordinates(x_positions, y_positions, z_positions):
    data = {
        "Координата X (м)": x_positions,
        "Координата Y (м)": y_positions,
        "Координата Z (м)": z_positions
    }
    df = pd.DataFrame(data)
    print("\nКоординаты спутника в системе ECEF (первые 10 точек):")
    print(df.head(10).to_string(index=False))  # вывод первых 10 строк
    return df


# основной код для обработки спутниковых данных
rinex_file_path = 'C:/Users/МАРТА/Desktop/Brdc0020.25n'  # путь к RINEX-файлу
satellite_id = 'G01'  # ID спутника, данные которого необходимо обработать
satellite_records = parse_rinex_file(rinex_file_path, satellite_id)  # чтение и извлечение данных
params = extract_orbital_params(satellite_records)  # извлечение параметров орбиты
x_positions, y_positions, z_positions = calculate_orbit(*params) # вычисление координат орбиты спутника
display_orbit_coordinates(x_positions, y_positions, z_positions) # вывод координат спутника в виде таблицы
plot_satellite_orbit(x_positions, y_positions, z_positions) # визуализация орбиты спутника