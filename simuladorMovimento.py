import numpy as np
import matplotlib.pyplot as plt

ace = 11  # Aceleração
vel = 40   # velocidade
posição_x = 30  
posição_y = 50
k = 3 # Resistencia
massa = 4    # Massa

def calcular_aceleracao(t, estado, k, m):
    x, y, vx, vy = estado
    v = np.sqrt(vx**2 + vy**2)
    ace_x = -k * v * vx
    ace_y = -k * v * vy
    ax = ace_x / m
    ay = (ace_y - m * ace) / m
    return [vx, vy, ax, ay]

def metodo_runge_kutta(t, estado, dt, k, m):
    k1 = calcular_aceleracao(t, estado, k, m)
    k2 = calcular_aceleracao(t + 0.5*dt, [s + 0.5*dt*k for s, k in zip(estado, k1)], k, m)
    k3 = calcular_aceleracao(t + 0.5*dt, [s + 0.5*dt*k for s, k in zip(estado, k2)], k, m)
    k4 = calcular_aceleracao(t + dt, [s + dt*k for s, k in zip(estado, k3)], k, m)
  
    return [s + dt/6 * (k1_i + 2*k2_i + 2*k3_i + k4_i) for s, k1_i, k2_i, k3_i, k4_i in zip(estado, k1, k2, k3, k4)]

def velocidade_com_resistencia(v0, theta, h0, k, m):
    if k==0 and theta==90:
        vx0=0
        vy0=v0
    else:
        theta = np.radians(theta)
        vx0 = v0 * np.cos(theta)
        vy0 = v0 * np.sin(theta)
    estado = [0, h0, vx0, vy0]
    dt = 0.01
    tempo = [0]
    x_array = [0]
    y_array = [h0]
    
    while estado[1] >= 0:
        t = tempo[-1]
        estado = metodo_runge_kutta(t, estado, dt, k, m)
        tempo.append(t + dt)
        x_array.append(estado[0])
        y_array.append(estado[1])
    
    tempo_aceleração = tempo[-1]
    
    return tempo_aceleração, estado, x_array, y_array


tempo_aceleração, tempo, x_array, y_array = \
velocidade_com_resistencia(vel, posição_x, posição_y, k, massa)

plt.plot(x_array, y_array)
plt.xlabel('Distância horizontal ')
plt.ylabel('Distância Vertical ')
plt.title('Resistencia')
plt.ylim(0,max(y_array)+7)
plt.plot([0,0],[0,posição_y],'k',linewidth=3)
plt.text(0.1,posição_y/2, "%g m" %(posição_y), fontsize=12, color='red')
plt.grid(True)
plt.show()

print(f'Tempo total : {tempo_aceleração:.2f}')