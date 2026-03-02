#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

struct Node {
    double x, y;
};

struct Element {
    int node_ids[4];
};

struct MacierzGlobalna {
    int size;
    double H_glob[16][16];

    MacierzGlobalna(int n) : size(n) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                H_glob[i][j] = 0.0;
            }
        }
    }

    void dodajLokalna(double H[4][4], const Element& element) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int global_i = element.node_ids[i] - 1;
                int global_j = element.node_ids[j] - 1;
                H_glob[global_i][global_j] += H[i][j];
            }
        }
    }

    void wypisz() const {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                cout << fixed << setprecision(5) << setw(10) << H_glob[i][j] << "\t";
            }
            cout << endl;
        }
    }
};

void wypiszMacierz(double macierz[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << fixed << setprecision(5) << setw(10) << macierz[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

struct RozwiazUkladRownan {

};

void obliczPochodne(double ksi, double eta, double dN_dKsi[], double dN_dEta[]) {
    dN_dKsi[0] = -0.25 * (1 - eta);
    dN_dKsi[1] = 0.25 * (1 - eta);
    dN_dKsi[2] = 0.25 * (1 + eta);
    dN_dKsi[3] = -0.25 * (1 + eta);

    dN_dEta[0] = -0.25 * (1 - ksi);
    dN_dEta[1] = -0.25 * (1 + ksi);
    dN_dEta[2] = 0.25 * (1 + ksi);
    dN_dEta[3] = 0.25 * (1 - ksi);
}

void obliczJakobian(double x[], double y[], double dN_dKsi[], double dN_dEta[], double J[2][2]) {
    J[0][0] = dN_dKsi[0] * x[0] + dN_dKsi[1] * x[1] + dN_dKsi[2] * x[2] + dN_dKsi[3] * x[3];
    J[0][1] = dN_dKsi[0] * y[0] + dN_dKsi[1] * y[1] + dN_dKsi[2] * y[2] + dN_dKsi[3] * y[3];
    J[1][0] = dN_dEta[0] * x[0] + dN_dEta[1] * x[1] + dN_dEta[2] * x[2] + dN_dEta[3] * x[3];
    J[1][1] = dN_dEta[0] * y[0] + dN_dEta[1] * y[1] + dN_dEta[2] * y[2] + dN_dEta[3] * y[3];
}

double wyznacznikJakobianu(double J[2][2]) {
    return J[0][0] * J[1][1] - J[0][1] * J[1][0];
}

void obliczH_pc(double H[4][4], double dN_dx[], double dN_dy[], double k) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            H[i][j] = k * (dN_dx[i] * dN_dx[j] + dN_dy[i] * dN_dy[j]);
        }
    }
}

void obliczHbc(double Hbc[4][4], const Node nodes[], const Element& element, double alpha) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Hbc[i][j] = 0.0;
        }
    }

    for (int edge = 0; edge < 4; edge++) {
        int node1 = element.node_ids[edge];
        int node2 = element.node_ids[(edge + 1) % 4];

        double x1 = nodes[node1 - 1].x;
        double y1 = nodes[node1 - 1].y;
        double x2 = nodes[node2 - 1].x;
        double y2 = nodes[node2 - 1].y;

        double length = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
        double detJ = length / 2.0;

        double N[2][4] = {
            {0.5, 0.5, 0.0, 0.0},
            {0.0, 0.5, 0.5, 0.0}
        };

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Hbc[i][j] += alpha * N[0][i] * N[0][j] * detJ;
            }
        }
    }
}

void obliczH(double x[], double y[], double H[4][4], double k, int liczbaPunktowCalkowania) {
    double punktyCalkowania[4], wagi[4];

    if (liczbaPunktowCalkowania == 2) {
        punktyCalkowania[0] = -1.0 / sqrt(3.0);
        punktyCalkowania[1] = 1.0 / sqrt(3.0);
        wagi[0] = wagi[1] = 1.0;
    }
    else if (liczbaPunktowCalkowania == 3) {
        punktyCalkowania[0] = -sqrt(3.0 / 5.0);
        punktyCalkowania[1] = 0.0;
        punktyCalkowania[2] = sqrt(3.0 / 5.0);
        wagi[0] = wagi[2] = 5.0 / 9.0;
        wagi[1] = 8.0 / 9.0;
    }
    else if (liczbaPunktowCalkowania == 4) {
        punktyCalkowania[0] = -sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0);
        punktyCalkowania[1] = -sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0);
        punktyCalkowania[2] = -punktyCalkowania[1];
        punktyCalkowania[3] = -punktyCalkowania[0];
        wagi[0] = wagi[3] = (18.0 - sqrt(30.0)) / 36.0;
        wagi[1] = wagi[2] = (18.0 + sqrt(30.0)) / 36.0;
    }
    else {
        cerr << "Niewspierana liczba punktow calkowania!" << endl;
        return;
    }

    for (int i = 0; i < liczbaPunktowCalkowania; i++) {
        for (int j = 0; j < liczbaPunktowCalkowania; j++) {
            double ksi = punktyCalkowania[i];
            double eta = punktyCalkowania[j];
            double waga = wagi[i] * wagi[j];

            double dN_dKsi[4], dN_dEta[4];
            obliczPochodne(ksi, eta, dN_dKsi, dN_dEta);

            double J[2][2];
            obliczJakobian(x, y, dN_dKsi, dN_dEta, J);
            double detJ = wyznacznikJakobianu(J);

            double J_odw[2][2] = {
                { J[1][1] / detJ, -J[0][1] / detJ },
                {-J[1][0] / detJ,  J[0][0] / detJ }
            };

            double dN_dx[4], dN_dy[4];
            for (int m = 0; m < 4; m++) {
                dN_dx[m] = J_odw[0][0] * dN_dKsi[m] + J_odw[0][1] * dN_dEta[m];
                dN_dy[m] = J_odw[1][0] * dN_dKsi[m] + J_odw[1][1] * dN_dEta[m];
            }

            double H_pc[4][4];
            obliczH_pc(H_pc, dN_dx, dN_dy, k * detJ * waga);

            for (int m = 0; m < 4; m++) {
                for (int n = 0; n < 4; n++) {
                    H[m][n] += H_pc[m][n];
                }
            }
        }
    }
}

bool wczytajWspolrzedne(Node nodes[], int liczbaWezlow, ifstream& plik) {
    string line;
    while (getline(plik, line)) {
        if (line.find("*Node") != string::npos) break;
    }
    for (int i = 0; i < liczbaWezlow; i++) {
        int id;
        char comma;
        plik >> id >> comma >> nodes[i].x >> comma >> nodes[i].y;
        if (plik.fail()) {
            cerr << "Błąd podczas wczytywania danych węzła o ID " << id << endl;
            return false;
        }
    }
    return true;
}

bool wczytajElementy(Element elements[], int liczbaElementow, ifstream& plik) {
    string line;
    while (getline(plik, line)) {
        if (line.find("*Element") != string::npos) break;
    }
    for (int i = 0; i < liczbaElementow; i++) {
        int id;
        char comma;
        plik >> id >> comma >> elements[i].node_ids[0] >> comma >> elements[i].node_ids[1]
            >> comma >> elements[i].node_ids[2] >> comma >> elements[i].node_ids[3];
        if (plik.fail()) {
            cerr << "Błąd podczas wczytywania danych elementu o ID " << id << endl;
            return false;
        }
    }
    return true;
}

int main() {
    const int liczbaWezlow = 16;
    const int liczbaElementow = 9;

    Node nodes[liczbaWezlow];
    Element elements[liczbaElementow];

    ifstream plik("Test2_4_4_MixGrid.txt");
    if (!plik) {
        cerr << "Nie udało się otworzyć pliku." << endl;
        return 1;
    }

    if (!wczytajWspolrzedne(nodes, liczbaWezlow, plik)) return 1;
    if (!wczytajElementy(elements, liczbaElementow, plik)) return 1;

    double alpha = 25.0;
    double k = 25.0;
    int liczbaPunktowCalkowania = 2;

    for (int e = 0; e < liczbaElementow; e++) {
        double x[4], y[4];
        for (int i = 0; i < 4; i++) {
            int node_index = elements[e].node_ids[i] - 1;
            x[i] = nodes[node_index].x;
            y[i] = nodes[node_index].y;
        }

        double H[4][4] = { 0 };
        obliczH(x, y, H, k, liczbaPunktowCalkowania);

        cout << "Macierz H dla elementu " << e + 1 << ":\n";
        wypiszMacierz(H);

        double Hbc[4][4] = { 0 };
        obliczHbc(Hbc, nodes, elements[e], alpha);

        cout << "Macierz HBC dla elementu " << e + 1 << ":\n";
        wypiszMacierz(Hbc);
    }

    return 0;
}
