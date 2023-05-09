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
    public partial class Form5 : Form
    {
        public Form5(double[,] ARR, int n, int m)
        {
            InitializeComponent();


            cPoint3D[,] i_Points3D = new cPoint3D[ARR.GetLength(0), ARR.GetLength(1)];
            for (int X = 0; X < ARR.GetLength(0); X++)
            {
                for (int Y = 0; Y < ARR.GetLength(1); Y++)
                {
                    i_Points3D[X, Y] = new cPoint3D(X * 2.0 / (double)n - 1, Y * 2.0 / (double)m - 1, ARR[X, Y]);
                }
            }

            // Setting one of the strings = null results in hiding this legend
            graph3D1.AxisX_Legend = "X from -1 to 1";
            graph3D1.AxisY_Legend = "Y from -1 to 1";
            graph3D1.AxisZ_Legend = "U-V";

            // IMPORTANT: Normalize X,Y,Z separately because there is an extreme mismatch
            // between X values (< 300) and Z values (> 30000)

            graph3D1.SetSurfacePoints(i_Points3D, eNormalize.Separate);

            graph3D1.Raster = (eRaster)3;
            Color[] c_Colors = Plot3D.ColorSchema.GetSchema((eSchema)12);
            graph3D1.SetColorScheme(c_Colors, 4);

        }
    }
}
