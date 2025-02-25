#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib> // Для вызова system

using namespace std;

// Функция для чтения данных из файла result.txt
vector<vector<double>> read_data(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }

    vector<vector<double>> data;
    string line;

    while (getline(file, line)) {
        vector<double> row;
        stringstream ss(line);
        double value;

        while (ss >> value) {
            row.push_back(value);
        }

        if (!row.empty()) {
            data.push_back(row);
        }
    }

    file.close();
    return data;
}

// Функция для создания временного файла данных для GNUPLOT
void create_temp_data_file(const vector<vector<double>>& data, const string& temp_filename) {
    ofstream temp_file(temp_filename);
    if (!temp_file.is_open()) {
        cerr << "Error: Cannot create temporary file " << temp_filename << endl;
        exit(1);
    }

    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[i].size(); ++j) {
            temp_file << i << " " << j << " " << data[i][j] << endl;
        }
        temp_file << endl; // Пустая строка между строками данных
    }

    temp_file.close();
}

int main() {
    // Имя файла с результатами
    string result_filename = "result.txt";
    string temp_data_filename = "temp_data.txt";

    // Чтение данных из файла
    vector<vector<double>> data = read_data(result_filename);

    // Создание временного файла данных для GNUPLOT
    create_temp_data_file(data, temp_data_filename);

    // Команда для GNUPLOT
    string gnuplot_command = R"(
        set pm3d map
        set palette rgbformulae 33,13,10
        set xrange [0:]
        set yrange [0:]
        splot ')" + temp_data_filename + R"(' using 1:2:3 with pm3d notitle
    )";

    // Запись команды в скрипт GNUPLOT
    string gnuplot_script = "gnuplot_script.tmp";
    ofstream script_file(gnuplot_script);
    if (!script_file.is_open()) {
        cerr << "Error: Cannot create GNUPLOT script file" << endl;
        exit(1);
    }

    script_file << gnuplot_command;
    script_file.close();

    // Вызов GNUPLOT
    cout << "Plotting data using GNUPLOT..." << endl;
    system(("gnuplot -persist " + gnuplot_script).c_str());

    // Удаление временных файлов
    remove(temp_data_filename.c_str());
    remove(gnuplot_script.c_str());

    return 0;
}