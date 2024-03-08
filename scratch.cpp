#include <vector>
using namespace std;
std::vector<std::vector<double>> convertToMatrixColumnWise(const std::vector<double>& vec, int n) {
    int N = vec.size();
    // Calculate the number of columns. If N is not perfectly divisible by n, add 1 more column to fit all elements.
    int m = (N % n == 0) ? (N / n) : (N / n + 1);

    std::vector<std::vector<double>> matrix(n, std::vector<double>(m, 0)); // Initialize matrix with zeros

    for (int i = 0; i < N; ++i) {
        // Calculate the row and column index based on column-wise filling
        int row = i % n;
        int col = i / n;
        matrix[row][col] = vec[i];
    }

    return matrix;
}



int main() {
    vector<double> vec = { 1,2,3,4,5,6,7,8,9,10 };
    vector<vector<double>> matrix = convertToMatrixColumnWise(vec, 3);
    matrix;

	return 0;
}