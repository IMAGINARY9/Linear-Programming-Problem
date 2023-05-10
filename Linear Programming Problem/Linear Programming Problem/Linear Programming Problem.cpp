#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

using namespace std;

typedef vector<vector<double>> Matrix;

void print_table(const Matrix& table) {
	int m = table.size(), n = table[0].size();

	// header
	cout << "+";
	for (int j = 0; j < n + 1; j++) {
		cout << "--------+";
	}
	cout << endl;
	cout << "|" << setw(3) << "i" << setw(6) << "|";
	for (int i = 0; i < n; i++)
	{
		cout << setw(3) << fixed << "X" << i << ((i > 9) ? setw(4) : setw(5)) << "|";
	}
	cout << endl;
	cout << "+";
	for (int j = 0; j < n + 1; j++) {
		cout << "--------+";
	}
	cout << endl;

	// table
	for (int i = 0; i < m; i++) {
		if (i + 1 == m) cout << "|" << setw(3) << "dj" << setw(6) << "|";
		else cout << "|" << setw(3) << i << setw(6) << "|";

		for (int j = 0; j < n; j++) {
			cout << setw(6) << fixed << setprecision(2) << table[i][j] << "  |";
		}
		cout << endl;
		cout << "+";
		for (int j = 0; j < n + 1; j++) {
			cout << "--------+";
		}
		cout << endl;
	}
	cout << endl;
}

void pivot(Matrix& table, int pivot_row, int pivot_col) {
	int m = table.size(), n = table[0].size();
	double pivot_val = table[pivot_row][pivot_col];
	cout << "pivot_val: " << pivot_val << endl << endl;

	for (int i = 1; i < m; i++) {
		if (i != pivot_row) {
			double mult = table[i][pivot_col];
			for (int j = 0; j < n; j++) {
				table[i][j] -= (mult * table[pivot_row][j]) / pivot_val;
			}
		}
	}

	for (int j = 0; j < n; j++) {
		table[pivot_row][j] /= pivot_val;
	}
}

int get_entering_col(const Matrix& table, bool max) {
	int last_row = table.size() - 1;
	int n = table[0].size();
	int entering_col = -1;
	double min_coeff = INFINITY;
	for (int j = 1; j < n - 1; j++) {
		double coeff = table[last_row][j];
		if (max)
		{
			if (coeff < min_coeff) {
				min_coeff = coeff;
				entering_col = j;
				if (min_coeff == 0) entering_col = -1;
			}
		}
		else {
			if (coeff < min_coeff && coeff > 0) {
				min_coeff = coeff;
				entering_col = j;
			}
		}
	}
	return entering_col;
}

int get_entering_row(const Matrix& table, bool max) {
	int last_row = table.size() - 1;
	int entering_row = -1;
	double min_coeff = INFINITY;
	for (int i = 1; i < last_row; i++)
	{
		double coeff = table[i][0];
		if (coeff < min_coeff && coeff != 0)
		{
			min_coeff = coeff;
			entering_row = i;
		}
	}

	return entering_row;
}

int get_leaving_col(const Matrix& table, int entering_row) {
	int last_row = table.size() - 1;
	int sz = table[0].size();
	int leaving_col = -1;
	double min_ratio = INFINITY;
	for (int i = 1; i < sz; i++) {
		double akj = table[entering_row][i];
		if (akj >= 0) continue;
		double ratio = -table[last_row][i] / akj;
		if (ratio < min_ratio) {
			min_ratio = ratio;
			leaving_col = i;
		}
	}

	return leaving_col;
}

int get_leaving_row(const Matrix& table, int entering_col) {
	int sz = table.size();
	int leaving_row = 1;
	double min_ratio = table[1][0] / table[1][entering_col];
	for (int i = 2; i < sz; i++) {
		double ratio = table[i][0] / table[i][entering_col];
		if (ratio < min_ratio && ratio > 0) {
			min_ratio = ratio;
			leaving_row = i;
		}
	}
	return leaving_row;
}

void results(Matrix& table, double c) {

	int m = table.size(), n = table[0].size();

	for (int i = 1; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (table[i][j] == 1)
			{
				for (int k = 0; k < m; k++)
					if (table[k][j] != 0) break;
				cout << 'x' << j << " = " << table[i][0] << endl;
			}
		}
	}

	cout << "F = " << table[m - 1][0] + c << endl;

}

Matrix first_it(Matrix& table, bool max) {

	int last_row = table.size() - 1;
	int sz = table[0].size();
	int r = 0;
	double min = INFINITY;

	for (int i = 0; i < sz; i++)
	{
		table[last_row][i] = table[last_row - 1][i] * 100 - table[0][i];
		min = (table[last_row][i] > 0 && table[last_row][i] < min)
			? table[last_row][i] : min;
		if (min == table[last_row][i]) r = i;
		table[last_row][i] = table[last_row - 1][i] * 1 - table[0][i];
	}

	int entering_col = r;
	if (max) entering_col = table[0].size() - 1;
	int leaving_row = get_leaving_row(table, entering_col);

	pivot(table, leaving_row, entering_col);
	cout << "1 iteration" << endl;
	print_table(table);

	return table;
}

Matrix simplex_method(Matrix& table, double c, bool max) {
	int counter = 2;
	while (true)
	{
		int entering_col = get_entering_col(table, max);


		if (entering_col == -1)
		{
			results(table, c);
			return table;
		}

		int leaving_row = get_leaving_row(table, entering_col);

		pivot(table, leaving_row, entering_col);
		cout << counter << " iteration" << endl;
		print_table(table);
		counter++;
	}
}

Matrix new_sol(Matrix& table)
{
	int last_row = table.size() - 1;
	int sz = table[0].size();
	vector<double> new_lim;
	int min_r = 0;
	double min_coeff = INFINITY;
	double max_num = 0;
	for (int i = 1; i < last_row; i++)
	{
		double coeff = table[i][0];
		coeff -= floor(coeff);
		if (coeff > max_num && coeff != 0)
		{
			max_num = coeff;
			min_r = i;
		}
	}

	for (int i = 0; i < last_row + 1; i++)
		table[i].push_back(0);

	for (int i = 0; i < sz; i++)
	{
		double num = table[min_r][i];
		num -= floor(num);
		new_lim.push_back(-num);
	}
	new_lim.push_back(1);

	table.push_back(new_lim);
	swap(table[last_row + 1], table[last_row]);
	print_table(table);
	return table;
}

bool checkX_round(Matrix& table)
{
	for (int i = 0; i < table.size(); i++)
	{
		double num = (round(table[i][0] * 100) / 100);
		if (num != round(num))
			return false;
	}
	return true;
}

bool checkX_positive(Matrix& table)
{
	for (int i = 0; i < table.size(); i++)
	{
		double num = (round(table[i][0] * 100) / 100);
		if (num < 0)
			return false;
	}
	return true;
}

Matrix dsimplex_method(Matrix& table, double c, bool max) {
	int counter = 1;
	while (true)
	{
		if (checkX_positive(table))
		{
			if (checkX_round(table))
			{
				results(table, c);
				return table;
			}
			table = new_sol(table);
		}

		int entering_row = get_entering_row(table, max);
		int leaving_col = get_leaving_col(table, entering_row);

		pivot(table, entering_row, leaving_col);
		cout << counter << " iteration" << endl;
		print_table(table);
		counter++;
	}
}

int main() {

	bool isGomoriMethod;
	cout << "You need to use Gomori method?\t1 \\ 0\n";
	cin >> isGomoriMethod;

	if (!isGomoriMethod)
	{
		bool isMax = false;

		Matrix table = {
			{0, -3, 1, 2, 0, 0, 0},
			{6, 1, 1, 2, 1, 0, 0},
			{3, -1, 1, -1, 0, -1, 1},
			{0, 0, 0, 0, 0, 0, 0}
		};
		double cnst = 11;

		print_table(table);

		auto matrix = first_it(table, isMax);
		simplex_method(matrix, cnst, isMax);

	}
	else
	{
		bool isMax = true;

		Matrix table = {
			{0, 8, 6, 0, 0},
			{12, 2, 5, 1, 0},
			{10, 4, 1, 0, 1},
			{0, 0, 0, 0, 0}
		};
		int cnst = 0;

		print_table(table);

		auto matrix = first_it(table, isMax);
		auto smatrix = simplex_method(matrix, cnst, isMax);
		auto new_sol_matrix = new_sol(smatrix);
		dsimplex_method(new_sol_matrix, cnst, isMax);
	}

	return 0;
}
