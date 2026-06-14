#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> 
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <sstream>
#include <string>
struct Particula {
    double m;                     // atribuição de massa unitária 
    double x, y, z;              // X Y Z
    double x_ant, y_ant, z_ant; // posição anterior
    double ax, ay, az;         // aceleração atual
};

void inicializar(Particula p[], int N, double dt) {
    for (int i = 0; i < N; i++) {
        // condicões de contorno 


p[i].x = (rand() / (double)RAND_MAX) - 0.5;
p[i].y = (rand() / (double)RAND_MAX) - 0.5;
p[i].z = (rand() / (double)RAND_MAX) - 0.5;
        p[i].m = 1.0;

        
        double vx = ((rand() / (double)RAND_MAX) - 0.5) * 0.01;
        double vy = 0; //momento linear 
        double vz = 0;

        
        p[i].x_ant = p[i].x - vx * dt;
        p[i].y_ant = p[i].y - vy * dt;
        p[i].z_ant = p[i].z - vz * dt;
    }
}



void calcularAceleracao(Particula p[], int N, double G, double epsilon) {
    for (int i = 0; i < N; i++) {
        p[i].ax = p[i].ay = p[i].az = 0; 
        for (int j = 0; j < N; j++) {
            if (i == j) continue;

            double dx = p[j].x - p[i].x;
            double dy = p[j].y - p[i].y;
            double dz = p[j].z - p[i].z;

            
            double r2 = dx*dx + dy*dy + dz*dz + epsilon*epsilon;
            double r3 = r2 * sqrt(r2); 

            
            double forca_comum = G * p[j].m / r3; // força por unidade de massa 
            p[i].ax += forca_comum * dx;
            p[i].ay += forca_comum * dy;
            p[i].az += forca_comum * dz;
        }
    }
}


void integrarVerlet(Particula p[], int N, double dt) {
    for (int i = 0; i < N; i++) {
        // Guarda a posição atual antes de atualizar
        double x_temp = p[i].x;
        double y_temp = p[i].y;
        double z_temp = p[i].z;

        //verlet evolução
        p[i].x = 2 * p[i].x - p[i].x_ant + p[i].ax * dt * dt;
        p[i].y = 2 * p[i].y - p[i].y_ant + p[i].ay * dt * dt;
        p[i].z = 2 * p[i].z - p[i].z_ant + p[i].az * dt * dt;

        // Atualiza o passado para o próximo passo
       
        p[i].x_ant = x_temp;
        p[i].y_ant = y_temp;
        p[i].z_ant = z_temp;
    }
}

double obterMassaTotal(Particula p[], int N) {
    double soma = 0;
    for (int i = 0; i < N; i++) {
        soma += p[i].m; 
    }
    return soma; //soma a massa total do sistema constante 
}


    void estabilizarVirial(Particula p[], int N, double G, double dt, double epsilon) {
    double K = 0.0; // Energia Cinética 
    double U = 0.0; // Energia Potencial 

    
    for (int i = 0; i < N; i++) {
        
        double vx = (p[i].x - p[i].x_ant) / dt;
        double vy = (p[i].y - p[i].y_ant) / dt;
        double vz = (p[i].z - p[i].z_ant) / dt;
        K += 0.5 * p[i].m * (vx*vx + vy*vy + vz*vz); //enrgia cinetica total
    }

    
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx = p[j].x - p[i].x;
            double dy = p[j].y - p[i].y;
            double dz = p[j].z - p[i].z;
           double r = sqrt(dx*dx + dy*dy + dz*dz + epsilon*epsilon); // softening
            U -= G * p[i].m * p[j].m / r;
        }
    }

    
    
    if (K > 0) {
        double fator = sqrt(fabs(U) / (2.0 * K));

        
        for (int i = 0; i < N; i++) {
            double vx = (p[i].x - p[i].x_ant) / dt * fator;
            double vy = (p[i].y - p[i].y_ant) / dt * fator;
            double vz = (p[i].z - p[i].z_ant) / dt * fator;

            p[i].x_ant = p[i].x - vx * dt;
            p[i].y_ant = p[i].y - vy * dt;
            p[i].z_ant = p[i].z - vz * dt;
        }
    }
}

double calcularTemperatura(Particula p[], int N, double dt) {
    double energiaCineticaTotal = 0;
    for (int i = 0; i < N; i++) {
        double vx = (p[i].x - p[i].x_ant) / dt;
        double vy = (p[i].y - p[i].y_ant) / dt;
        double vz = (p[i].z - p[i].z_ant) / dt;
        energiaCineticaTotal += 0.5 * p[i].m * (vx*vx + vy*vy + vz*vz); //energia micorscopica 
    }
    
    return (2.0 * energiaCineticaTotal) / (3.0 * N); //media por particula (temperatura estastistica)
}
// Para o cálculo do raio de meia-massa
struct Diagnostico {
    double energiaTotal;
    double energiaCinetica;
    double energiaPotencial;
    double temperatura;
    double raioMeiaMassa;
};

        Diagnostico calcularTermodinamica(Particula p[], int N, double G, double dt, double epsilon) {
    double K = 0.0; // Energia Cinética
    double U = 0.0; // Energia Potencial
    double cmX = 0, cmY = 0, cmZ = 0;
    double massaTotal = 0;

    
    for (int i = 0; i < N; i++) {
        massaTotal += p[i].m;
        cmX += p[i].x * p[i].m; 
        cmY += p[i].y * p[i].m;
        cmZ += p[i].z * p[i].m;
                    //referencia do sisteema 
        double vx = (p[i].x - p[i].x_ant) / dt;
        double vy = (p[i].y - p[i].y_ant) / dt;
        double vz = (p[i].z - p[i].z_ant) / dt;
        K += 0.5 * p[i].m * (vx*vx + vy*vy + vz*vz);
    }
    cmX /= massaTotal; cmY /= massaTotal; cmZ /= massaTotal;

    
    std::vector<double> raios;
    for (int i = 0; i < N; i++) {
       
        double r_cm = sqrt(pow(p[i].x - cmX, 2) + pow(p[i].y - cmY, 2) + pow(p[i].z - cmZ, 2));
        raios.push_back(r_cm);

        for (int j = i + 1; j < N; j++) {
            double dx = p[j].x - p[i].x;
            double dy = p[j].y - p[i].y;
            double dz = p[j].z - p[i].z;
            double r = sqrt(dx*dx + dy*dy + dz*dz + epsilon*epsilon); // Softening
            U -= G * p[i].m * p[j].m / r; //energia potecial total
        }
    }

    
    std::sort(raios.begin(), raios.end());
    double rhalf = raios[N / 2];
    double T = (2.0 * K) / (3.0 * N);
    return { 
    (K + U),  // energia total
    K,        // energia cinética
    U,        // energia potencial
    T,        // temperatura
    rhalf     // raio
};
}


void carregarEstado(Particula p[], int &N, double &dt, std::ifstream &arquivo) {
    std::string linha;

    
    while (std::getline(arquivo, linha)) {
        if (linha == "# estado_inicial") break;
    }

    
    std::getline(arquivo, linha);

    int i = 0;
    while (std::getline(arquivo, linha)) {
        if (linha.empty() || linha[0] == '#') break;

        std::stringstream ss(linha);
        std::string valor;

        std::getline(ss, valor, ','); // id

        std::getline(ss, valor, ',');
        p[i].x = std::stod(valor);

        std::getline(ss, valor, ',');
        p[i].y = std::stod(valor);

        std::getline(ss, valor, ',');
        p[i].z = std::stod(valor);

        double vx, vy, vz;

        std::getline(ss, valor, ',');
        vx = std::stod(valor);

        std::getline(ss, valor, ',');
        vy = std::stod(valor);

        std::getline(ss, valor, ',');
        vz = std::stod(valor);

        std::getline(ss, valor, ',');
        p[i].m = std::stod(valor);

        
        p[i].x_ant = p[i].x - vx * dt;
        p[i].y_ant = p[i].y - vy * dt;
        p[i].z_ant = p[i].z - vz * dt;

        i++;
    }

    N = i;
}


int main() {
    const int N = 10;
    const double G = 1.0;
    const double epsilon = 0.01;
    const double dt = 0.01;
    
    Particula p1[N];
    Particula p2[N];

srand(123);
inicializar(p1, N, dt);

// copia p1 --Ç> p2
for (int i = 0; i < N; i++) {
    p2[i] = p1[i];
}

// pequena perturbação  para criar caos  
p2[0].x += 1e-6;

estabilizarVirial(p1, N, G, dt, epsilon);
    estabilizarVirial(p2, N, G, dt, epsilon);

    std::ofstream arquivo("caos.csv");
    
    arquivo << "step,tempo,distancia,E,K,U,T,raio,virial\n";
    

    //loop
    for (int step = 0; step < 100000; step++) {

        double tempo = step * dt;

        calcularAceleracao(p1, N, G, epsilon);
        integrarVerlet(p1, N, dt);

        calcularAceleracao(p2, N, G, epsilon);
        integrarVerlet(p2, N, dt);

        if (step % 100 == 0) {

        Diagnostico d1 = calcularTermodinamica(p1, N, G, dt, epsilon);
            double virial = 2*d1.energiaCinetica + d1.energiaPotencial;
            double distancia = 0;

    for (int i = 0; i < N; i++) {
        double dx = p1[i].x - p2[i].x;
        double dy = p1[i].y - p2[i].y;
        double dz = p1[i].z - p2[i].z;
        distancia += dx*dx + dy*dy + dz*dz;
    }

    distancia = sqrt(distancia / N);


        arquivo << step << ","
        << tempo << ","
        << distancia << ","
        << d1.energiaTotal << ","
        << d1.energiaCinetica << ","
        << d1.energiaPotencial << ","
        << d1.temperatura << ","
        << d1.raioMeiaMassa << ","
        << virial
        << "\n";
        }
    }
     arquivo.close(); 
    
    return 0;
}