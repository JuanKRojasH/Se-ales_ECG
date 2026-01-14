"""

Analizador Ritmo Cardiaco

Código completamente reestructurado y optimizado para análisis ECG

"""

import numpy as np                  # Operaciones numéricas y manejo de arrays
import matplotlib.pyplot as plt     # Generación de gráficos y visualizaciones
import scipy.io as sciio            # Carga y manejo de archivos .MAT
import scipy.signal as sig          # Procesamiento de señales digitales

# ============================================================================
# SECCIÓN 1: CARGA Y PREPARACIÓN DE DATOS CARDIACOS
# ============================================================================

archivo_principal = "datos_cardiologicos_brutos.mat"       # Definir nombre del archivo principal
contenido_principal = sciio.loadmat(archivo_principal)     # Cargar datos desde archivo .MAT
array_senal_ecg = np.double(contenido_principal['senal'])  # Extraer y convertir señal ECG
frecuencia = np.double(contenido_principal['frecuencia'])  # Extraer frecuencia de muestreo
array_senal_ecg = array_senal_ecg[0]                       # Aplanar array a una dimensión

archivo_referencia = "plantilla_qrs_estandar.mat"          # Definir archivo de patrón referencia
contenido_referencia = sciio.loadmat(archivo_referencia)   # Cargar patrón de referencia
plantilla_qrs = np.double(contenido_referencia['plantilla']) # Extraer plantilla QRS
frecuencia_ref = np.double(contenido_referencia['frecuencia']) # Extraer frecuencia referencia
plantilla_qrs = plantilla_qrs[0]                           # Aplanar array de plantilla

# ============================================================================
# SECCIÓN 2: CONSTRUCCIÓN DE ESCALA TEMPORAL
# ============================================================================

numero_muestras = np.size(array_senal_ecg)                 # Calcular total de muestras disponibles
escala_temporal = np.arange(0.0, numero_muestras) / frecuencia # Crear vector tiempo en segundos

# ============================================================================
# SECCIÓN 3: VISUALIZACIÓN INICIAL DE DATOS
# ============================================================================

plt.figure(figsize=(14, 4))                                # Crear figura para visualización
plt.plot(escala_temporal, array_senal_ecg, 'b-', alpha=0.7, linewidth=0.8) # Graficar señal ECG
plt.xlabel("Tiempo (segundos)", fontsize=11)               # Etiqueta eje X
plt.ylabel("Amplitud (mV)", fontsize=11)                   # Etiqueta eje Y
plt.title("Registro Electrocardiográfico Original", fontsize=13) # Título del gráfico
plt.grid(True, alpha=0.3)                                  # Activar cuadrícula transparente
plt.tight_layout()                                         # Ajustar márgenes automáticamente

plt.figure(figsize=(10, 3))                                # Nueva figura para patrón
plt.plot(plantilla_qrs, 'r-', linewidth=1.5)               # Graficar plantilla de referencia
plt.xlabel("Índice de Muestra", fontsize=11)               # Etiqueta eje X
plt.ylabel("Amplitud Normalizada", fontsize=11)            # Etiqueta eje Y
plt.title("Plantilla de Complejo QRS de Referencia", fontsize=13) # Título
plt.grid(True, alpha=0.3)                                  # Cuadrícula transparente
plt.tight_layout()                                         # Ajustar márgenes

# ============================================================================
# SECCIÓN 4: PROCESAMIENTO CON OPERADORES DIFERENCIALES
# ============================================================================

senal_diferenciada = np.diff(array_senal_ecg)              # Calcular primera derivada señal ECG
plantilla_diferenciada = np.diff(plantilla_qrs)            # Calcular primera derivada plantilla
plantilla_invertida = np.flip(plantilla_diferenciada)      # Invertir plantilla en dominio temporal

# ============================================================================
# SECCIÓN 5: CÁLCULO DE CORRELACIÓN CRUZADA
# ============================================================================

correlacion_resultante = np.convolve(senal_diferenciada, plantilla_invertida, mode='same') # Correlación cruzada
correlacion_umbralizada = correlacion_resultante * (correlacion_resultante > 0.0)          # Umbralizar valores positivos

# ============================================================================
# SECCIÓN 6: DETECCIÓN AUTOMÁTICA DE EVENTOS CARDIACOS
# ============================================================================

indices_picos_brutos = sig.find_peaks(correlacion_umbralizada)[0]        # Detectar picos iniciales
valor_maximo = np.max(correlacion_umbralizada[indices_picos_brutos])     # Calcular máximo global
umbral_inferior = valor_maximo * 0.25                                    # Definir umbral mínimo adaptativo
separacion_minima = int(np.round(0.25 * frecuencia))                     # Calcular separación mínima (250ms)

indices_picos_filtrados = sig.find_peaks(correlacion_umbralizada,        # Detección refinada con criterios
                                         height=[umbral_inferior, valor_maximo],
                                         distance=separacion_minima)[0]

# ============================================================================
# SECCIÓN 7: CÁLCULO DE INTERVALOS RR Y FRECUENCIA CARDIACA
# ============================================================================

intervalos_rr = np.diff(indices_picos_filtrados) / frecuencia            # Calcular intervalos RR en segundos
frecuencia_cardiaca_instantanea = 60.0 / intervalos_rr                   # Convertir a latidos por minuto
frecuencia_promedio = np.mean(frecuencia_cardiaca_instantanea)           # Calcular frecuencia promedio
variabilidad_rr = np.std(intervalos_rr) * 1000                           # Calcular variabilidad en milisegundos

# ============================================================================
# SECCIÓN 8: VISUALIZACIÓN COMPARATIVA DE RESULTADOS
# ============================================================================

fig, axes = plt.subplots(3, 1, figsize=(15, 9), sharex=True)            # Crear figura con 3 subgráficos

axes[0].plot(escala_temporal[1:], senal_diferenciada, 'teal', linewidth=0.8, alpha=0.8) # Señal diferenciada
axes[0].set_ylabel("ECG Derivado", fontsize=11)                         # Etiqueta eje Y primer panel
axes[0].grid(True, alpha=0.2)                                           # Cuadrícula sutil
axes[0].set_title("Análisis Completo de Detección de Complejos QRS", fontsize=14) # Título principal

axes[1].plot(escala_temporal[1:], correlacion_umbralizada, 'forestgreen', linewidth=1.0) # Señal correlación
axes[1].plot(escala_temporal[indices_picos_filtrados], 
            correlacion_umbralizada[indices_picos_filtrados], 
            'ro', markersize=6, alpha=0.7)                              # Marcar picos detectados
axes[1].set_ylabel("Correlación", fontsize=11)                          # Etiqueta eje Y segundo panel
axes[1].grid(True, alpha=0.2)                                           # Cuadrícula sutil

axes[2].plot(escala_temporal[indices_picos_filtrados[1:]], 
            frecuencia_cardiaca_instantanea, 
            'darkorange', linewidth=1.5, marker='o', markersize=4)      # Graficar frecuencia cardíaca
axes[2].set_xlabel("Tiempo (segundos)", fontsize=11)                    # Etiqueta eje X
axes[2].set_ylabel("Frecuencia (lpm)", fontsize=11)                     # Etiqueta eje Y tercer panel
axes[2].axhline(y=frecuencia_promedio, color='r', linestyle='--', alpha=0.5, label=f'Promedio: {frecuencia_promedio:.1f}') # Línea promedio
axes[2].legend(loc='upper right', fontsize=10)                          # Leyenda
axes[2].grid(True, alpha=0.2)                                           # Cuadrícula sutil
axes[2].set_xlim([escala_temporal[0], escala_temporal[1200]])           # Limitar vista a primeros 1200 puntos

plt.tight_layout()                                                      # Ajuste automático de márgenes

# ============================================================================
# SECCIÓN 9: REPORTE ESTADÍSTICO Y EXPORTACIÓN
# ============================================================================

print("\n" + "="*65)                                                    # Línea separadora
print("RESUMEN ESTADÍSTICO DEL ANÁLISIS CARDIACO")                     # Encabezado
print("="*65)                                                          # Línea separadora
print(f"• Total de latidos detectados: {len(indices_picos_filtrados)}") # Conteo latidos
print(f"• Frecuencia cardíaca promedio: {frecuencia_promedio:.1f} lpm") # Frecuencia promedio
print(f"• Frecuencia cardíaca mínima: {np.min(frecuencia_cardiaca_instantanea):.1f} lpm") # Mínimo
print(f"• Frecuencia cardíaca máxima: {np.max(frecuencia_cardiaca_instantanea):.1f} lpm") # Máximo
print(f"• Variabilidad RR (SDNN): {variabilidad_rr:.1f} ms")           # Variabilidad
print(f"• Duración análisis: {escala_temporal[-1]:.1f} segundos")       # Duración total
print("="*65)                                                          # Línea separadora

# Exportar resultados a archivos de texto
np.savetxt('indices_latidos_detectados.txt', indices_picos_filtrados, fmt='%d', header='Índices de latidos detectados') # Guardar índices
np.savetxt('intervalos_rr_calculados.txt', intervalos_rr, fmt='%.4f', header='Intervalos RR en segundos') # Guardar intervalos
np.savetxt('frecuencia_cardiaca_instantanea.txt', frecuencia_cardiaca_instantanea, fmt='%.2f', header='Frecuencia cardiaca instantánea (lpm)') # Guardar frecuencias

print("\nArchivos exportados exitosamente:")                           # Mensaje confirmación
print("  1. indices_latidos_detectados.txt")                          # Lista archivo 1
print("  2. intervalos_rr_calculados.txt")                            # Lista archivo 2
print("  3. frecuencia_cardiaca_instantanea.txt")                     # Lista archivo 3
print("\nProceso completado satisfactoriamente ✓")                     # Mensaje final

# Generar histograma de frecuencia cardíaca (análisis adicional)
plt.figure(figsize=(10, 5))                                           # Nueva figura para histograma
plt.hist(frecuencia_cardiaca_instantanea, bins=15, color='steelblue', edgecolor='black', alpha=0.7) # Histograma
plt.xlabel("Frecuencia Cardíaca (latidos por minuto)", fontsize=11)   # Etiqueta eje X
plt.ylabel("Frecuencia de Ocurrencia", fontsize=11)                   # Etiqueta eje Y
plt.title("Distribución de Frecuencia Cardíaca Instantánea", fontsize=13) # Título
plt.grid(True, alpha=0.3)                                             # Cuadrícula
plt.axvline(x=frecuencia_promedio, color='red', linestyle='--', linewidth=2, label=f'Media: {frecuencia_promedio:.1f} lpm') # Línea media
plt.legend()                                                          # Mostrar leyenda
plt.tight_layout()                                                    # Ajustar márgenes
plt.show()                                                            # Mostrar gráfico final
