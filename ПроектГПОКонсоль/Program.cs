using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ПроектГПОКонсоль
{
    delegate double TFunc(double[] x);

    class Program
    {
        static void Main(string[] args)
        {
            int SearchAgents_no, Max_iteration, dim; string Function_name;
            double[] lb; double[] ub; TFunc fobj; double lb1, ub1;
            SearchAgents_no = 40; // Number of search agents
            Function_name = "F1"; // Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
            Max_iteration = 500; // Maximum numbef of iterations
            //Load details of the selected benchmark function
            Get_Functions_details(Function_name, out lb1, out ub1, out dim, out fobj);
            DA(SearchAgents_no, Max_iteration, lb, ub, dim, fobj);
        }
        static void Get_Functions_details(string F, out double lb1, out double ub1, out int dim, out TFunc fobj)
        {
            fobj = F1;
            lb1 = -100;
            ub1 = 100;
            dim = 10; 
            switch (F)
            {
                case "F1":
                    break;
                case "F2":
                    fobj = F2;
                    lb1 = -10;
                    ub1 = 10;
                    dim = 10;
                    break;
                case "F3":
                    fobj = F3;
                    lb1 = -100;
                    ub1 = 100;
                    dim = 10;
                    break;
                case "F4":
                    fobj = F4;
                    lb1 = -100;
                    ub1 = 100;
                    dim = 10;
                    break;
                case "F5":
                    fobj = F5;
                    lb1 = -30;
                    ub1 = 30;
                    dim = 10;
                    break;
                case "F6":
                    fobj = F6;
                    lb1 = -100;
                    ub1 = 100;
                    dim = 10;
                    break;
                case "F7":
                    fobj = F7;
                    lb1 = -1.28;
                    ub1 = 1.28;
                    dim = 10;
                    break;
                case "F8":
                    fobj = F8;
                    lb1 = -500;
                    ub1 = 500;
                    dim = 10;
                    break;
                case "F9":
                    fobj = F9;
                    lb1 = -5.12;
                    ub1 = 5.12;
                    dim = 10;
                    break;
                case "F10":
                    fobj = F10;
                    lb1 = -32;
                    ub1 = 32;
                    dim = 10;
                    break;
                case "F11":
                    fobj = F11;
                    lb1 = -600;
                    ub1 = 600;
                    dim = 10;
                    break;
                case "F12":
                    fobj = F12;
                    lb1 = -50;
                    ub1 = 50;
                    dim = 10;
                    break;
                case "F13":
                    fobj = F13;
                    lb1 = -50;
                    ub1 = 50;
                    dim = 10;
                    break;
            }
        }
        static double F1(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + Math.Pow(elem, 2);
            }
            return sum;
        }
        static double F2(double[] x)
        {
            double sum = 0, prod = 1; foreach (double elem in x)
            {
                sum = sum + Math.Abs(elem);
                prod = prod * Math.Abs(elem);
            }
            return (sum + prod);
        }
        static double F3(double[] x)
        {
            int i, j, dim = x.Length; double y = 0;
            for (i = 0; i < dim; ++i)
            {
                double sum = 0;
                for (j = 0; j <= i; ++j)
                {
                    sum = sum + x[j];
                }
                y = y + Math.Pow(sum, 2);
            }
            return y;
        }
        static double F4(double[] x)
        {
            double max = 0;
            foreach (double elem in x)
            {
                if (max < Math.Abs(elem)) max = Math.Abs(elem);
            }
            return max;
        }
        static double F5(double[] x)
        {
            int dim = x.Length, i; double sum = 0;
            for (i = 0; i <= (dim - 2); ++i)
            {
                sum = sum + (100 * Math.Pow((x[i + 1] - Math.Pow(x[i], 2)), 2) + Math.Pow((x[i] - 1), 2));
            }
            return sum;
        }
        static double F6(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + Math.Pow((elem + 0.5), 2);
            }
            return sum;
        }
        static double F7(double[] x)
        {
            Random rnd = new Random();
            double sum = rnd.NextDouble(); int i, dim = x.Length;
            for (i = 0; i < dim; ++i)
            {
                sum = sum + (i + 1) * Math.Pow(x[i], 4);
            }
            return sum;
        }
        static double F8(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + (-1 * elem * Math.Sin(Math.Sqrt(Math.Abs(elem))));
            }
            return sum;
        }
        static double F9(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + (Math.Pow(elem, 2) - 10 * Math.Cos(2 * Math.PI * elem) + 10);
            }
            return sum;
        }
        static double F10(double[] x)
        {
            double sum1 = 0, sum2 = 0, sum;
            foreach (double elem in x)
            {
                sum1 = sum1 + Math.Pow(elem, 2);
                sum2 = sum2 + Math.Cos(2 * Math.PI * elem);
            }
            sum1 = sum1 / x.Length; sum2 = sum2 / x.Length;
            sum = -20 * Math.Exp(-0.2 * Math.Sqrt(sum1)) - Math.Exp(sum2) + 20 + Math.Exp(1);
            return sum;
        }
        static double F11(double[] x)
        {
            double sum = 0, mply = 1; int i, dim = x.Length;
            for (i = 0; i < dim; ++i)
            {
                sum = sum + Math.Pow(x[i], 2);
                mply = mply * Math.Cos(x[i] / Math.Sqrt(i + 1));
            }
            double y = sum / 4000 - mply + 1;
            return y;
        }
        static double F12(double[] x)
        {
            int dim = x.Length;
            double[] z = new double[dim];
            for (int i = 0; i < dim; ++i)
            {
                z[i] = 1 + (x[i] + 1) / 4.0;
            }
            double sum1 = 0, sum2 = 0;
            for (int i = 0; i < (dim - 1); ++i)
            {
                sum1 = sum1 + Math.Pow(z[i] - 1, 2) * (1 + 10 * Math.Pow(Math.Sin(Math.PI * z[i + 1]), 2));
            }
            foreach (double elem in x)
            {
                sum2 = sum2 + uFunc(elem, 10, 100, 4);
            }
            double y = (Math.PI / dim) * (10 * Math.Sin(Math.PI * z[0]) + sum1 + Math.Pow(z[dim - 1] - 1, 2)) + sum2;
            return y;
        }
        static double uFunc(double x, int a, int k, int m)
        {
            if (x > a) return (k * Math.Pow(x - a, m));
            if (x < -a) return (k * Math.Pow(-x - a, m));
            return 0;
        }
        static double F13(double[] x)
        {
            int dim = x.Length; double sum1 = 0, sum2 = 0;
            foreach (double elem in x)
            {
                sum1 = sum1 + Math.Pow(elem - 1, 2) * (1 + Math.Pow(Math.Sin(3 * Math.PI * elem + 1), 2));
                sum2 = sum2 + uFunc(elem, 5, 100, 4);
            }
            double y = 0.1 * (Math.Pow(Math.Sin(3 * Math.PI * x[0]), 2) + sum1 + Math.Pow(x[dim - 1] - 1, 2) * (1 + Math.Pow(Math.Sin(2 * Math.PI * x[dim - 1]), 2))) + sum2;
            return y;
        }
        static void initialization(int SearchAgents_no, int dim, double[] lb, double[] ub, double[,] X)
        {
            int i, j; Random rnd = new Random();
            for (i = 0; i < dim; ++i)
            {
                for (j = 0; j < SearchAgents_no; ++j)
                {
                    X[i, j] = rnd.NextDouble() * (ub[i] - lb[i]) + lb[i];
                }
            }
        }
        static void DA(int SearchAgents_no, int Max_iteration, double[] lb, double[] ub, int dim, TFunc fobj)
        {

        }
    }
}