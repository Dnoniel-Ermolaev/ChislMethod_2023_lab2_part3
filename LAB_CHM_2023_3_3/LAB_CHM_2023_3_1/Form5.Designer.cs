
namespace LAB_CHM_2023_3_1
{
    partial class Form5
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.graph3D1 = new Plot3D.Graph3D();
            this.SuspendLayout();
            // 
            // graph3D1
            // 
            this.graph3D1.AxisX_Color = System.Drawing.Color.DarkBlue;
            this.graph3D1.AxisX_Legend = null;
            this.graph3D1.AxisY_Color = System.Drawing.Color.DarkGreen;
            this.graph3D1.AxisY_Legend = null;
            this.graph3D1.AxisZ_Color = System.Drawing.Color.DarkRed;
            this.graph3D1.AxisZ_Legend = null;
            this.graph3D1.BackColor = System.Drawing.Color.White;
            this.graph3D1.BorderColor = System.Drawing.Color.FromArgb(((int)(((byte)(180)))), ((int)(((byte)(180)))), ((int)(((byte)(180)))));
            this.graph3D1.Location = new System.Drawing.Point(13, 13);
            this.graph3D1.Name = "graph3D1";
            this.graph3D1.PolygonLineColor = System.Drawing.Color.Black;
            this.graph3D1.Raster = Plot3D.Graph3D.eRaster.Off;
            this.graph3D1.Size = new System.Drawing.Size(810, 599);
            this.graph3D1.TabIndex = 0;
            this.graph3D1.TopLegendColor = System.Drawing.Color.FromArgb(((int)(((byte)(200)))), ((int)(((byte)(200)))), ((int)(((byte)(150)))));
            // 
            // Form5
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(835, 624);
            this.Controls.Add(this.graph3D1);
            this.Name = "Form5";
            this.Text = "МОДУЛЬ РАЗНОСТИ ЗНАЧЕНИЙ ГРАФИКОВ";
            this.ResumeLayout(false);

        }

        #endregion

        private Plot3D.Graph3D graph3D1;
    }
}