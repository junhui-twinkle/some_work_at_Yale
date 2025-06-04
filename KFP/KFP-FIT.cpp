#undef UNICODE
#undef _UNICODE
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<graphics.h>
#include<conio.h>

#define max_points 100
#define error_threshold1 0.0001
#define error_threshold2 0.0001

FILE* open_file(const char* filename);
int read_data(FILE* file, double* t, double* ratio, const char* filename);
double model1(double alpha, double k_x, double t);
double model2(double alpha, double beta, double k_x, double k_y, double t);
double error_function1(double alpha, double k_x, double* t, double* ratio, int number_count);
double error_function2(double alpha, double beta, double k_x, double k_y, double* t, double* ratio, int number_count);
double gradient_descent1(double* alpha, double* k_x, double* t, double* ratio, int number_count, double learning_rate, int iterations);
double gradient_descent2(double* alpha, double* beta, double* k_x, double* k_y, double* t, double* ratio, int number_count, double learning_rate, int iterations);
double R_square1(double alpha, double k_x, double* t, double* ratio, int number_count);
double R_square2(double alpha, double beta, double k_x, double k_y, double* t, double* ratio, int number_count);
double R_square1_adj(double alpha, double k_x, double* t, double* ratio, int number_count);
double R_square2_adj(double alpha, double beta, double k_x, double k_y, double* t, double* ratio, int number_count);
void coordinate_axis(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double* t, int number_count);
void plot_fitted_curve1(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double alpha, double k_x, double* t, double* ratio, int number_count);
void plot_fitted_curve2(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double alpha, double beta, double k_x, double k_y, double* t, double* ratio, int number_count);
void plot_data(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double* t, double* ratio, int number_count);
void figure1(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double alpha, double k_x, double* t, double* ratio, int number_count);
void figure2(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double alpha, double beta, double k_x, double k_y, double* t, double* ratio, int number_count);
errno_t output_file(const char* filename, double alpha, double beta, double k_x, double k_y, double Rsquare1, double Rsquare2, double Rsquare1_adj, double Rsquare2_adj);

int main()
{
	FILE* file1, *file2;
	double t1[max_points], t2[max_points], ratio1[max_points], ratio2[max_points];
	int number_count1 = 0, number_count2 = 0, max_iterations1 = 10000, max_iterations2 = 20000, output_status = 0;
	int width = 930, height = 670, width_edge = 80, height_edge = 60, n_x = 11, n_y = 11;
	double alpha = 1, k_x = 1, beta = 1, k_y = 1, learning_rate1 = 1, learning_rate2 = 1, gradient_descent11 = 0,gradient_descent22 = 0;
	double Rsquare1 = 0, Rsquare2 = 0, Rsquare1_adj = 0, Rsquare2_adj = 0;
	char user_input[100];
	char output_filename[100];
	printf("Please enter the output filename (without extension).\n");
	scanf_s("%s", user_input, 100);
	sprintf_s(output_filename, "%s.txt", user_input);
	file1 = open_file("d1.txt");
	if (file1 == NULL)
	{
		return EXIT_FAILURE;
	}
	file2 = open_file("d2.txt");
	if (file2 == NULL)
	{
		return EXIT_FAILURE;
	}
	number_count1 = read_data(file1, t1, ratio1, "d1.txt");
	if (number_count1 == -1)
	{
		return EXIT_FAILURE;
	}
	number_count2 = read_data(file2, t2, ratio2, "d2.txt");
	if (number_count2 == -1)
	{
		return EXIT_FAILURE;
	}
	fclose(file1);
	fclose(file2);
	gradient_descent11 = gradient_descent1(&alpha, &k_x, t1, ratio1, number_count1, learning_rate1, max_iterations1);
	gradient_descent22 = gradient_descent2(&alpha, &beta, &k_x, &k_y, t2, ratio2, number_count2, learning_rate2, max_iterations2);
	printf("alpha = %lf\nk_x = %lf\n", alpha, k_x);
	printf("beta = %lf\nk_y = %lf\n", beta, k_y);
	Rsquare1 = R_square1(alpha, k_x, t1, ratio1, number_count1);
	Rsquare2 = R_square2(alpha, beta, k_x, k_y, t2, ratio2, number_count2);
	Rsquare1_adj = R_square1_adj(alpha, k_x, t1, ratio1, number_count1);
	Rsquare2_adj = R_square2_adj(alpha, beta, k_x, k_y, t2, ratio2, number_count2);
	printf("R1 = %lf\nR2 = %lf\n", Rsquare1, Rsquare2);
	printf("R1_adj = %lf\nR2_adj = %lf\n", Rsquare1_adj, Rsquare2_adj);
	figure1(width, height, width_edge, height_edge, n_x, n_y, alpha, k_x, t1, ratio1, number_count1);
	figure2(width, height, width_edge, height_edge, n_x, n_y, alpha, beta, k_x, k_y, t2, ratio2, number_count2);
	output_status = output_file(output_filename, alpha, beta, k_x, k_y, Rsquare1, Rsquare2, Rsquare1_adj, Rsquare2_adj);
	if (output_status == 0)
	{
		printf("Parameters writtern to file %s successfully.\n", output_filename);
	}
	else
		printf("Failed to write to the file %s.\n", output_filename);
	return EXIT_SUCCESS;
}

FILE* open_file(const char* filename)
{
	FILE* file;
	errno_t err = fopen_s(&file, filename, "r");
	if (err != 0)
	{
		perror("Error opening file");
		return NULL;
	}
	return file;
}

//the fourth parameter is used to output more informative messages
int read_data(FILE* file, double* t, double* ratio, const char* filename)
{
	int number_count = 0;
	while (1)
	{
		if (fscanf_s(file, "%lf %lf", &t[number_count], &ratio[number_count]) == 2 && number_count < max_points)
		{
			number_count++;
		}
		else
		{
			break;
		}
	}
	if (feof(file))
	{
		printf("Successfully read %d points from the file %s.\n", number_count, filename);
	}
	else
	{
		perror("Error reading file1");
		fclose(file);
		return EOF;
	}
	return number_count;
}

double model1(double alpha, double k_x, double t)
{
	return (1 - alpha) * exp(-k_x * t) + alpha;
}

double model2(double alpha, double beta, double k_x, double k_y, double t)
{
	return ((1 - alpha) * (1 - beta) / (k_x - k_y)) * (k_x * exp(-k_y * t) - k_y * exp(-k_x * t)) + (1 - (1 - alpha) * (1 - beta));
}

double error_function1(double alpha, double k_x, double* t, double* ratio, int number_count)
{
	double error = 0;
	for (int i = 0; i < number_count; i++)
	{
		error += pow(ratio[i] - model1(alpha, k_x, t[i]), 2);
	}
	return error / number_count;
}

double error_function2(double alpha, double beta, double k_x, double k_y, double* t, double* ratio, int number_count)
{
	double error = 0;
	for (int i = 0; i < number_count; i++)
	{
		error += pow(ratio[i] - model2(alpha, beta, k_x, k_y, t[i]), 2);
	}
	return error / number_count;
}

double gradient_descent1(double* alpha, double* k_x, double* t, double* ratio, int number_count, double learning_rate, int max_iterations)
{
	double previous_err = error_function1(*alpha, *k_x, t, ratio, number_count);
	double current_err = previous_err;
	int count = 0;
	for (int i = 0; i < max_iterations; i++)
	{
		double grad_alpha = 0;
		double grad_k_x = 0;
		for (int j = 0; j < number_count; j++)
		{
			double predict = model1(*alpha, *k_x, t[j]);
			double difference = predict - ratio[j];
			grad_alpha += 2 * difference * (-exp(-(*k_x) * t[j]) + 1);
			grad_k_x += 2 * difference * (*alpha - 1) * exp(-(*k_x) * t[j]) * t[j];
		}
		*alpha -= learning_rate * grad_alpha / number_count;
		*k_x -= learning_rate * grad_k_x / number_count;
		if (isnan(*alpha) || isnan(*k_x))
		{
			printf("Encountered NaN during training at iteration %d.\n", i + 1);
			break;
		}
		current_err = error_function1(*alpha, *k_x, t, ratio, number_count);
		//printf("iteration: %d previous error = %lf\tcurrent error %lf\t difference = %lf\n", i + 1, previous_err, current_err, previous_err - current_err);
		if (fabs(current_err) < error_threshold1)
		{
			printf("Converge.\n");
			break;
		}
		if (count >= 10)
		{
			printf("Stable.\n");
			break;
		}
		if (current_err == previous_err)
		{
			count++;
		}
		else
		{
			count = 0;
		}
		previous_err = current_err;
	}
	return current_err;
}

double gradient_descent2(double* alpha, double* beta, double* k_x, double* k_y, double* t, double* ratio, int number_count, double learning_rate, int max_iterations)
{
	double previous_err = error_function2(*alpha, *beta, *k_x, *k_y, t, ratio, number_count);
	double current_err = previous_err;
	int count = 0;
	double epsilon = 0.00001;
	for (int i = 0; i < max_iterations; i++)
	{
		double grad_beta = 0;
		double grad_k_y = 0;
		for (int j = 0; j < number_count; j++)
		{
			double predict = model2(*alpha, *beta, *k_x, *k_y, t[j]);
			double difference = predict - ratio[j];
			/*double nominator = *k_x - (*k_y);
			if (fabs(nominator) < 0.06)
			{
				nominator = 0.06;
			}*/
			double err_beta_plus = error_function2(*alpha, *beta + epsilon, *k_x, *k_y, t, ratio, number_count);
			double err_beta_miuns = error_function2(*alpha, *beta - epsilon, *k_x, *k_y, t, ratio, number_count);
			double err_k_y_plus = error_function2(*alpha, *beta, *k_x, *k_y + epsilon, t, ratio, number_count);
			double err_k_y_miuns = error_function2(*alpha, *beta, *k_x, *k_y - epsilon, t, ratio, number_count);
			grad_beta = (err_beta_plus - err_beta_miuns) / (2 * epsilon);
			grad_k_y = (err_k_y_plus - err_k_y_miuns) / (2 * epsilon);
			//grad_beta += 2 * difference * ((*k_x * exp(-(*k_y) * t[j]) - *k_y * exp(-(*k_x) * t[j])) * (*alpha - 1) / nominator + 1 - *alpha);
			//grad_k_y += 2 * difference * (1 - (*alpha)) * (1 - (*beta)) / pow(nominator, 2) * (*k_x * exp(-(*k_y) * t[j]) - *k_y * exp(-(*k_x) * t[j])) + (1 - (*alpha)) * (1 - (*beta)) / (nominator) * (*k_x * (-t[j]) * exp(-(*k_y) * t[j]) - exp(-(*k_x) * t[j]));
		}
		*beta -= learning_rate * grad_beta / number_count;
		*k_y -= learning_rate * grad_k_y / number_count;
		if (isnan(*beta) || isnan(*k_y))
		{
			printf("Encountered NaN during training at iteration %d.\n", i + 1);
			break;
		}
		current_err = error_function2(*alpha, *beta, *k_x, *k_y, t, ratio, number_count);
		//printf("iteration: %d previous error = %lf\tcurrent error %lf\t difference = %lf\n", i + 1, previous_err, current_err, previous_err - current_err);
		//printf("alpha = %lf\t beta = %lf\t k_x = %lf\t k_y = %lf\n", *alpha, *beta, *k_x, *k_y);
		if (fabs(current_err) < error_threshold1)
		{
			printf("Converge.\n");
			break;
		}
		if (count >= 10)
		{
			printf("Stable.\n");
			break;
		}
		if (fabs(current_err - previous_err) < 1e-17)
		{
			count++;
		}
		else
		{
			count = 0;
		}
		previous_err = current_err;
	}
	return current_err;
}

double R_square1(double alpha, double k_x, double* t, double* ratio, int number_count)
{
	double sum = 0;
	double ratio_average = 0, SST =0, SSR =0;
	for (int i = 0; i < number_count; i++)
	{
		sum += ratio[i];
	}
	ratio_average = sum / number_count;
	for (int i = 0; i < number_count; i++)
	{
		SST += pow((ratio[i] - ratio_average), 2);
		SSR += pow((ratio[i] - model1(alpha, k_x, t[i])), 2);
	}
	return 1 - SSR / SST;
}

double R_square1_adj(double alpha, double k_x, double* t, double* ratio, int number_count)
{
	double Rsquare = R_square1(alpha, k_x, t, ratio, number_count);
	int p = 1;
	double Rsquare_adj = 1 - (1 - Rsquare) * (number_count - 1) / (number_count - p - 1);
	return Rsquare_adj;
}

double R_square2(double alpha, double beta, double k_x, double k_y, double* t, double* ratio, int number_count)
{
	double sum = 0;
	double ratio_average = 0, SST = 0, SSR = 0;
	for (int i = 0; i < number_count; i++)
	{
		sum += ratio[i];
	}
	ratio_average = sum / number_count;
	for (int i = 0; i < number_count; i++)
	{
		SST += pow((ratio[i] - ratio_average), 2);
		SSR += pow((ratio[i] - model2(alpha, beta, k_x, k_y, t[i])), 2);
	}
	return 1 - SSR / SST;
}

double R_square2_adj(double alpha, double beta, double k_x, double k_y, double* t, double* ratio, int number_count)
{
	double Rsquare = R_square2(alpha, beta, k_x, k_y, t, ratio, number_count);
	int p = 1;
	double Rsquare_adj = 1 - (1 - Rsquare) * (number_count - 1) / (number_count - p - 1);
	return Rsquare_adj;
}

void coordinate_axis(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double* t, int number_count)
{
	initgraph(width, height, 0);
	setbkcolor(WHITE);
	cleardevice();
	int x_tick_length = (getwidth() - 2 * width_edge) / n_x;
	int y_tick_length = (getheight() - 2 * height_edge) / n_y;
	double ratio_tick[10] = { 0, 0.20, 0.40, 0.60, 0.80, 1.00 };
	double t_tick[20] = { 0 };
	for (int i = 1; i < 6; i++)
	{
		t_tick[i] = i * t[number_count - 1] / 5;
	}
	setlinecolor(BLACK);
	line(width_edge, getheight() - height_edge, getwidth() - width_edge, getheight() - height_edge); // bottom
	line(width_edge, height_edge, getwidth() - width_edge, height_edge);//top
	line(width_edge, height_edge, width_edge, getheight() - height_edge);//left
	line(getwidth() - width_edge, height_edge, getwidth() - width_edge, getheight() - height_edge);//right
	settextcolor(BLACK);
	for (int i = 0; i < n_y + 1; i++)
	{
		line(width_edge - 10, getheight() - height_edge - i * y_tick_length, width_edge, getheight() - height_edge - i * y_tick_length);
		char ratio_tick_char[20];
		if (i % 2 == 0)
		{
			sprintf_s(ratio_tick_char, "%.2lf", ratio_tick[i / 2]);
			outtextxy(width_edge - 30, getheight() - height_edge - i * y_tick_length - 5, ratio_tick_char);
		}
	}
	for (int i = 0; i < n_x + 1; i++)
	{
		line(width_edge + i * x_tick_length, getheight() - height_edge, width_edge + i * x_tick_length, getheight() - height_edge + 10);
		char t_tick_char[20];
		if (i % 2 == 0)
		{
			sprintf_s(t_tick_char, "%.2lf", t_tick[i / 2]);
			outtextxy(width_edge + i * x_tick_length - 10, getheight() - height_edge + 20, t_tick_char);
		}
	}
}

void plot_fitted_curve1(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double alpha, double k_x, double* t, double* ratio, int number_count)
{
	setlinecolor(BLUE);
	setfillcolor(BLUE);
	int x_length = (getwidth() - 2 * width_edge);
	int y_length = (getheight() - 2 * height_edge);
	int x_start = width_edge;
	int x_end = width_edge + (n_x - 1) * x_length / n_x;
	int y_start = getheight() - height_edge;
	int y_end = getheight() - height_edge - (n_y - 1) * y_length / n_y;
	for (int i = 0; i < x_length * (n_x - 1) / n_x; i++)
	{
		int coordinate_x1 = width_edge + i;
		double t_x1 = (double)t[number_count - 1] / (x_end - x_start) * (coordinate_x1 - x_start);
		double ratio_y1 = model1(alpha, k_x, t_x1);
		int coordinate_y1 = (int)((y_end - y_start) * ratio_y1 + y_start);
		int coordinate_x2 = width_edge + i + 1;
		double t_x2 = (double)t[number_count - 1] / (x_end - x_start) * (coordinate_x2 - x_start);
		double ratio_y2 = model1(alpha, k_x, t_x2);
		int coordinate_y2 = (int)((y_end - y_start) * ratio_y2 + y_start);
		line(coordinate_x1, coordinate_y1, coordinate_x2, coordinate_y2);
	}
}

void plot_fitted_curve2(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double alpha, double beta, double k_x, double k_y, double* t, double* ratio, int number_count)
{
	setlinecolor(BLUE);
	int x_length = (getwidth() - 2 * width_edge);
	int y_length = (getheight() - 2 * height_edge);
	int x_start = width_edge;
	int x_end = width_edge + (n_x - 1) * x_length / n_x;
	int y_start = getheight() - height_edge;
	int y_end = getheight() - height_edge - (n_y - 1) * y_length / n_y;
	for (int i = 0; i < x_length * (n_x - 1) / n_x; i++)
	{
		int coordinate_x1 = width_edge + i;
		double t_x1 = (double)t[number_count - 1] / (x_end - x_start) * (coordinate_x1 - x_start);
		double ratio_y1 = model2(alpha, beta, k_x, k_y, t_x1);
		int coordinate_y1 = (int)((y_end - y_start) * ratio_y1 + y_start);
		int coordinate_x2 = width_edge + i + 1;
		double t_x2 = (double)t[number_count - 1] / (x_end - x_start) * (coordinate_x2 - x_start);
		double ratio_y2 = model2(alpha, beta, k_x, k_y, t_x2);
		int coordinate_y2 = (int)((y_end - y_start) * ratio_y2 + y_start);
		line(coordinate_x1, coordinate_y1, coordinate_x2, coordinate_y2);
	}
}

void plot_data(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double* t, double* ratio, int number_count)
{
	setlinecolor(RED);
	setfillcolor(RED);
	int radius = 3;
	int x_length = (getwidth() - 2 * width_edge);
	int y_length = (getheight() - 2 * height_edge);
	int x_start = width_edge;
	int x_end = width_edge + (n_x - 1) * x_length / n_x;
	int y_start = getheight() - height_edge;
	int y_end = getheight() - height_edge - (n_y - 1) * y_length / n_y;
	for (int i = 0; i < number_count; i++)
	{
		int coordinate_x1 = (int)((x_end - x_start) * t[i] / t[number_count - 1] + x_start);
		int coordinate_y1 = (int)((y_end - y_start) * ratio[i] + y_start);
		fillcircle(coordinate_x1, coordinate_y1, radius);
	}
}

void figure1(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double alpha, double k_x, double* t, double* ratio, int number_count)
{
	coordinate_axis(width, height, width_edge, height_edge, n_x, n_y, t, number_count);
	plot_fitted_curve1(width, height, width_edge, height_edge, n_x, n_y, alpha, k_x, t, ratio, number_count);
	plot_data(width, height, width_edge, height_edge, n_x, n_y, t, ratio, number_count);
	saveimage("figure1.bmp", NULL);
	_getch();
	closegraph();
}

void figure2(int width, int height, int width_edge, int height_edge, int n_x, int n_y, double alpha, double beta, double k_x, double k_y, double* t, double* ratio, int number_count)
{
	coordinate_axis(width, height, width_edge, height_edge, n_x, n_y, t, number_count);
	plot_fitted_curve2(width, height, width_edge, height_edge, n_x, n_y, alpha, beta, k_x, k_y, t, ratio, number_count);
	plot_data(width, height, width_edge, height_edge, n_x, n_y, t, ratio, number_count);
	saveimage("figure2.bmp", NULL);
	_getch();
	closegraph();
}

errno_t output_file(const char* filename, double alpha, double beta, double k_x, double k_y, double Rsquare1, double Rsquare2, double Rsquare1_adj, double Rsquare2_adj)
{
	FILE* file;
	errno_t err = fopen_s(&file, filename, "w");
	if (err != 0)
	{
		perror("Error opening file");
		return err;
	}
	fprintf(file, "alpha\t%lf\n", alpha);
	fprintf(file, "kx\t%lf\n", k_x);
	fprintf(file, "Rsquare\t%lf\n", Rsquare1);
	fprintf(file, "Rsquare_adj\t%lf\n", Rsquare1_adj);
	fprintf(file, "\n");
	fprintf(file, "beta\t%lf\n", beta);
	fprintf(file, "ky\t%lf\n", k_y);
	fprintf(file, "Rsquare\t%lf\n", Rsquare2);
	fprintf(file, "Rsquare_adj\t%lf\n", Rsquare2_adj);
	fclose(file);
	return 0;
}