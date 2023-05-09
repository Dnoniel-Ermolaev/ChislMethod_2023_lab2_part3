using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

using delRendererFunction = Plot3D.Graph3D.delRendererFunction;
using cPoint3D = Plot3D.Graph3D.cPoint3D;
using eRaster = Plot3D.Graph3D.eRaster;
using cScatter = Plot3D.Graph3D.cScatter;
using eNormalize = Plot3D.Graph3D.eNormalize;
using eSchema = Plot3D.ColorSchema.eSchema;

namespace LAB_CHM_2023_3_1
{
    public partial class Form4 : Form
    {
        double[,] ARR;
        int N;
        int M;
        public Form4(double[][] u, double[][] v1, int n, int m, bool r)
        {
            InitializeComponent();
          

            if (r == false)
            {
                ARR = new double[n + 1, m + 1];
                double[,] arr = new double[n + 1, m + 1];
                for (int i = 0; i < n+1; i++)
                {
                    for (int j = 0; j < m + 1; j++)
                    {
                        arr[i, j] = u[i][j];
                    }
                };

                double[,] arr2 = new double[n + 1, m + 1];
                for (int i = 0; i < n + 1; i++)
                {
                    for (int j = 0; j < m + 1; j++)
                    {
                        arr2[i, j] = v1[i][j];
                    }
                };

                cPoint3D[,] i_Points3D = new cPoint3D[arr.GetLength(0), arr.GetLength(1)];
                for (int X = 0; X < arr.GetLength(0); X++)
                {
                    for (int Y = 0; Y < arr.GetLength(1); Y++)
                    {
                        i_Points3D[X, Y] = new cPoint3D(X * 2.0 / (double)n - 1, Y * 2.0 / (double)m - 1, arr[X, Y]);
                    }
                }

                cPoint3D[,] i_Points3D2 = new cPoint3D[arr2.GetLength(0), arr2.GetLength(1)];
                for (int X = 0; X < arr2.GetLength(0); X++)
                {
                    for (int Y = 0; Y < arr2.GetLength(1); Y++)
                    {
                        i_Points3D2[X, Y] = new cPoint3D(X * 2.0 / (double)n - 1, Y * 2.0 / (double)m - 1, arr2[X, Y]);
                    }
                }
                // Setting one of the strings = null results in hiding this legend
                graph3D1.AxisX_Legend = "X from -1 to 1";
                graph3D1.AxisY_Legend = "Y from -1 to 1";
                graph3D1.AxisZ_Legend = "U";

                graph3D2.AxisX_Legend = "X from -1 to 1";
                graph3D2.AxisY_Legend = "Y from -1 to 1";
                graph3D2.AxisZ_Legend = "V";
                // IMPORTANT: Normalize X,Y,Z separately because there is an extreme mismatch
                // between X values (< 300) and Z values (> 30000)

                graph3D1.SetSurfacePoints(i_Points3D, eNormalize.Separate);

                graph3D2.SetSurfacePoints(i_Points3D2, eNormalize.Separate);

                graph3D1.Raster = (eRaster)3;
                Color[] c_Colors = Plot3D.ColorSchema.GetSchema((eSchema)0);
                graph3D1.SetColorScheme(c_Colors, 4);

                graph3D2.Raster = (eRaster)3;
                Color[] c_Colors2 = Plot3D.ColorSchema.GetSchema((eSchema)2);
                graph3D2.SetColorScheme(c_Colors2, 1);


                //ДЛЯ БУДУЩЕЙ РАЗНИЦЫ

                for (int i = 0; i < n + 1; i++)
                {
                    for (int j = 0; j < m + 1; j++)
                    {
                        ARR[i, j] = Math.Abs(arr[i, j] - arr2[i, j]);
                    }
                };
                N = n;
                M = m;

            }
            else
            {
                ARR = new double[n + 1, m + 1];

                double[,] arr = new double[n + 1, m + 1];
                for (int i = 0; i < n+1; i++)
                {
                    for (int j = 0; j < m + 1; j++)
                    {
                        arr[i, j] = u[i][j];
                    }
                };

                double[,] arr2 = new double[n*2 + 1, m*2 + 1];
                for (int i = 0; i < n*2 + 1; i++)
                {
                    for (int j = 0; j < m*2 + 1; j++)
                    {
                        arr2[i, j] = v1[i][j];
                    }
                };

                cPoint3D[,] i_Points3D = new cPoint3D[arr.GetLength(0), arr.GetLength(1)];
                for (int X = 0; X < arr.GetLength(0); X++)
                {
                    for (int Y = 0; Y < arr.GetLength(1); Y++)
                    {
                        i_Points3D[X, Y] = new cPoint3D(X * 2.0 / (double)n - 1, Y * 2.0 / (double)m - 1, arr[X, Y]);
                    }
                }

                cPoint3D[,] i_Points3D2 = new cPoint3D[arr2.GetLength(0), arr2.GetLength(1)];
                for (int X = 0; X < arr2.GetLength(0); X++)
                {
                    for (int Y = 0; Y < arr2.GetLength(1); Y++)
                    {
                        i_Points3D2[X, Y] = new cPoint3D(X * 2.0 / (double)n - 1, Y * 2.0 / (double)m - 1, arr2[X, Y]);
                    }
                }
                graph3D1.AxisX_Legend = "X from -1 to 1";
                graph3D1.AxisY_Legend = "Y from -1 to 1";
                graph3D1.AxisZ_Legend = "V";

                graph3D2.AxisX_Legend = "X from -1 to 1";
                graph3D2.AxisY_Legend = "Y from -1 to 1";
                graph3D2.AxisZ_Legend = "V_2";

                // IMPORTANT: Normalize X,Y,Z separately because there is an extreme mismatch
                // between X values (< 300) and Z values (> 30000)

                graph3D1.SetSurfacePoints(i_Points3D, eNormalize.Separate);

                graph3D2.SetSurfacePoints(i_Points3D2, eNormalize.Separate);

                graph3D1.Raster = (eRaster)3;
                Color[] c_Colors = Plot3D.ColorSchema.GetSchema((eSchema)0);
                graph3D1.SetColorScheme(c_Colors, 4);

                graph3D2.Raster = (eRaster)3;
                Color[] c_Colors2 = Plot3D.ColorSchema.GetSchema((eSchema)2);
                graph3D2.SetColorScheme(c_Colors2, 1);


                //ДЛЯ БУДУЩЕЙ РАЗНИЦЫ

                for (int i = 0; i < n + 1; i++)
                {
                    for (int j = 0; j < m + 1; j++)
                    {
                        ARR[i, j] = Math.Abs(arr[i, j] - arr2[i*2, j*2]);
                    }
                };
                N = n;
                M = m;
            }

        }

        private void SetFormula(string s_Formula)
        {
            try
            {
                delRendererFunction f_Function = Plot3D.FunctionCompiler.Compile(s_Formula);

                // IMPORTANT: Normalize maintainig the relation between X,Y,Z values otherwise the function will be distorted.
                graph3D1.SetFunction(f_Function, new PointF(-10, -10), new PointF(10, 10), 0.5, eNormalize.MaintainXYZ);
            }
            catch (Exception Ex)
            {
                MessageBox.Show(Ex.Message, "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }
        private void SetSurface(double[,] s32_Values)
        {
            cPoint3D[,] i_Points3D = new cPoint3D[s32_Values.GetLength(0), s32_Values.GetLength(1)];
            for (int X = 0; X < s32_Values.GetLength(0); X++)
            {
                for (int Y = 0; Y < s32_Values.GetLength(1); Y++)
                {
                    i_Points3D[X, Y] = new cPoint3D(X * 10, Y * 500, s32_Values[X, Y]);
                }
            }

            // IMPORTANT: Normalize X,Y,Z separately because there is an extreme mismatch
            // between X values (< 300) and Z values (> 30000)
            graph3D1.SetSurfacePoints(i_Points3D, eNormalize.Separate);
        }

        private void button1_Click(object sender, EventArgs e)
        {
            Form5 form5 = new Form5(ARR, N, M);
            form5.Show();
        }
    }
}
