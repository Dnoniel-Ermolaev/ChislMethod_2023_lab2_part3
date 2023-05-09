
namespace LAB_CHM_2023_3_1
{
    partial class Form4
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
            this.graph3D2 = new Plot3D.Graph3D();
            this.graph3D1 = new Plot3D.Graph3D();
            this.button1 = new System.Windows.Forms.Button();
            this.label1 = new System.Windows.Forms.Label();
            this.SuspendLayout();
            // 
            // graph3D2
            // 
            this.graph3D2.AxisX_Color = System.Drawing.Color.DarkBlue;
            this.graph3D2.AxisX_Legend = null;
            this.graph3D2.AxisY_Color = System.Drawing.Color.DarkGreen;
            this.graph3D2.AxisY_Legend = null;
            this.graph3D2.AxisZ_Color = System.Drawing.Color.DarkRed;
            this.graph3D2.AxisZ_Legend = null;
            this.graph3D2.BackColor = System.Drawing.Color.White;
            this.graph3D2.BorderColor = System.Drawing.Color.FromArgb(((int)(((byte)(180)))), ((int)(((byte)(180)))), ((int)(((byte)(180)))));
            this.graph3D2.Location = new System.Drawing.Point(512, 12);
            this.graph3D2.Name = "graph3D2";
            this.graph3D2.PolygonLineColor = System.Drawing.Color.Black;
            this.graph3D2.Raster = Plot3D.Graph3D.eRaster.Off;
            this.graph3D2.Size = new System.Drawing.Size(460, 504);
            this.graph3D2.TabIndex = 1;
            this.graph3D2.TopLegendColor = System.Drawing.Color.FromArgb(((int)(((byte)(200)))), ((int)(((byte)(200)))), ((int)(((byte)(150)))));
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
            this.graph3D1.Location = new System.Drawing.Point(12, 12);
            this.graph3D1.Name = "graph3D1";
            this.graph3D1.PolygonLineColor = System.Drawing.Color.Black;
            this.graph3D1.Raster = Plot3D.Graph3D.eRaster.Off;
            this.graph3D1.Size = new System.Drawing.Size(460, 504);
            this.graph3D1.TabIndex = 0;
            this.graph3D1.TopLegendColor = System.Drawing.Color.FromArgb(((int)(((byte)(200)))), ((int)(((byte)(200)))), ((int)(((byte)(150)))));
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(492, 522);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(75, 23);
            this.button1.TabIndex = 2;
            this.button1.Text = "РАСЧЁТ";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(370, 527);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(116, 13);
            this.label1.TabIndex = 3;
            this.label1.Text = "ВЫВЕСТИ РАЗНИЦУ";
            // 
            // Form4
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(984, 561);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.graph3D2);
            this.Controls.Add(this.graph3D1);
            this.Name = "Form4";
            this.Text = "ДВА ГРАФИКА ПО ТАБЛИЦАМ";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private Plot3D.Graph3D graph3D1;
        private Plot3D.Graph3D graph3D2;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Label label1;
    }
}