# Simulación de Energía Solar con Paneles Fotovoltaicos

## Descripción del proyecto
Este proyecto permite calcular la **posición solar**, la **irradiancia horizontal global (GHI)** y la **potencia generada por un panel fotovoltaico** a lo largo de un día completo usando Python.  

Se incluyen visualizaciones de la **altura solar** y la **potencia del panel** para cada hora del día.

El proyecto es útil para:
- Dimensionamiento de sistemas fotovoltaicos.
- Evaluación de la generación energética en distintas ubicaciones.
- Análisis educativo de geometría solar y radiación solar.

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta, timezone

# =====================================
# FUNCIONES AUXILIARES


# Función que devuelve el día del año (1-365/366)
def dia_del_año(dt): 
    return dt.timetuple().tm_yday

# Convierte la fecha y hora en una fracción del año en radianes
def fraccion_del_año(dt):
    n = dia_del_año(dt)
    hora = dt.hour + dt.minute/60 + dt.second/3600
    # Centrado en mediodía para mejorar precisión
    return 2*math.pi/365 * (n - 1 + (hora - 12)/24)

# Calcula la declinación solar (δ) según el día del año
def declinacion(gamma):
    return (0.006918 - 0.399912*math.cos(gamma) + 0.070257*math.sin(gamma)
            - 0.006758*math.cos(2*gamma) + 0.000907*math.sin(2*gamma)
            - 0.002697*math.cos(3*gamma) + 0.00148*math.sin(3*gamma))

# Ecuación del tiempo (minutos) para corregir la hora solar
def ecuacion_del_tiempo(gamma):
    return 229.18*(0.000075 + 0.001868*math.cos(gamma) - 0.032077*math.sin(gamma)
                   - 0.014615*math.cos(2*gamma) - 0.040849*math.sin(2*gamma))

# Calcula el ángulo horario H (°) según hora local, longitud y zona horaria
def angulo_horario(t, lon, utc):
    gamma = fraccion_del_año(t)
    EoT = ecuacion_del_tiempo(gamma)  # corrección de ecuación del tiempo
    LSTM = 15*utc  # hora solar estándar de la zona
    TC = EoT + 4*(LSTM - lon)   # corrección de tiempo total (minutos)
    minutos_locales = t.hour*60 + t.minute + t.second/60
    hora_solar = (minutos_locales + TC)/60
    return 15*(hora_solar - 12)  # Ángulo horario en grados

# Calcula altura solar (alt), azimut (az) y ángulo cenital (Z)
def posicion_solar(t, lat, lon, utc):
    phi = math.radians(lat)  # latitud en radianes
    delta = declinacion(fraccion_del_año(t))  # declinación solar
    H = math.radians(angulo_horario(t, lon, utc))  # ángulo horario
    # coseno del ángulo cenital
    cosZ = math.sin(phi)*math.sin(delta) + math.cos(phi)*math.cos(delta)*math.cos(H)
    cosZ = max(-1, min(1, cosZ))  # limitar a [-1,1] para acos
    Z = math.degrees(math.acos(cosZ))  # ángulo cenital en grados
    alt = 90 - Z  # altura solar sobre el horizonte
    alt_rad = math.radians(alt)
    # cálculo del azimut solar
    sinA = -math.sin(H)*math.cos(delta)/max(math.cos(alt_rad), 1e-6)
    cosA = (math.sin(delta) - math.sin(phi)*math.sin(alt_rad))/(math.cos(phi)*max(math.cos(alt_rad), 1e-6))
    az = (math.degrees(math.atan2(sinA, cosA)) + 360) % 360
    return alt, az, Z

# Modelo de irradiancia Haurwitz (W/m²) simplificado
def haurwitz(cosZ):
    return 0 if cosZ <= 0 else 1098*cosZ*math.exp(-0.057/max(cosZ,1e-6))

# =====================================
# ENTRADAS DEL USUARIO

lat = float(input("Latitud (°): "))
lon = float(input("Longitud (°): "))
fecha = input("Fecha (YYYY-MM-DD): ")
utc = int(input("Zona horaria UTC (ej: -5): "))
tilt = float(input("Inclinación del panel (°): "))
azim = float(input("Azimut del panel (°): "))
area = float(input("Área total del panel (m²): "))
eff = float(input("Eficiencia del panel (0-1): "))

# =====================================
# SIMULACIÓN HORARIA

# Crear rango de 24 horas para la fecha indicada
date_obj = datetime.strptime(fecha, "%Y-%m-%d")
tz = timezone(timedelta(hours=utc))
horas = pd.date_range(start=date_obj, periods=24, freq="H", tz=tz)

resultados = []

# Bucle de cálculo hora a hora
for t in horas:
    # 1. Posición solar
    alt, azs, Z = posicion_solar(t, lat, lon, utc)
    
    # 2. Irradiancia horizontal (GHI)
    cosZ = math.cos(math.radians(Z))
    GHI = haurwitz(cosZ)
    
    # 3. Ángulo de incidencia en el panel
    alt_r, azs_r, tilt_r, azim_r = map(math.radians, [alt, azs, tilt, azim])
    cosAOI = (math.sin(alt_r)*math.cos(tilt_r) +
              math.cos(alt_r)*math.sin(tilt_r)*math.cos(azs_r-azim_r))
    cosAOI = max(cosAOI, 0)
    
    # 4. Radiación sobre el panel (POA)
    POA = GHI * cosAOI
    
    # 5. Potencia generada por el panel (kW)
    Potencia = POA * area * eff / 1000
    
    # Guardar resultados
    resultados.append([t, alt, azs, GHI, POA, Potencia])

# Crear DataFrame con resultados
df = pd.DataFrame(resultados, columns=["Tiempo","Altura°","Azimut°","GHI(W/m²)","POA(W/m²)","Potencia(kW)"])
print(df)

# Energía diaria total
energia = df["Potencia(kW)"].sum()
print(f"\nEnergía diaria estimada: {energia:.2f} kWh")

# =====================================
# VISUALIZACIÓN DE RESULTADOS


# Altura solar vs hora
plt.figure(figsize=(9,5))
plt.plot(df["Tiempo"], df["Altura°"], marker="o", color="orange")
plt.fill_between(df["Tiempo"], df["Altura°"], color="orange", alpha=0.2)
plt.title("Altura solar vs Hora del día")
plt.xlabel("Hora")
plt.ylabel("Altura solar (°)")
plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=2))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
plt.grid(alpha=0.4)
plt.gcf().autofmt_xdate()
plt.savefig("altitud_solar.png", dpi=200)
plt.show()

# Potencia del panel vs hora
plt.figure(figsize=(9,5))
plt.plot(df["Tiempo"], df["Potencia(kW)"], marker="o", color="green", label="Potencia FV")
plt.fill_between(df["Tiempo"], df["Potencia(kW)"], color="green", alpha=0.2)
plt.title(f"Potencia del panel FV vs Hora del día\nEnergía total: {energia:.2f} kWh")
plt.xlabel("Hora")
plt.ylabel("Potencia (kW)")
plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=2))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
plt.grid(alpha=0.4)
plt.gcf().autofmt_xdate()
plt.legend()
plt.savefig("potencia_pv.png", dpi=200)
plt.show()

