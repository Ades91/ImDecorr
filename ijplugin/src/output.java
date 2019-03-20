/* 
% 
%   Output functions for the ImageDecorrelationAnalysis plugin
%
% ---------------------------------------
%
%   Copyright © 2018 Adrien Descloux - adrien.descloux@epfl.ch, 
%   École Polytechnique Fédérale de Lausanne, LBEN/LOB,
%   BM 5.134, Station 17, 1015 Lausanne, Switzerland.
%
%  	This program is free software: you can redistribute it and/or modify
%  	it under the terms of the GNU General Public License as published by
% 	the Free Software Foundation, either version 3 of the License, or
%  	(at your option) any later version.
%
%  	This program is distributed in the hope that it will be useful,
%  	but WITHOUT ANY WARRANTY; without even the implied warranty of
%  	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  	GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%  	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

import java.awt.Color;
import java.util.Calendar;

import ij.ImagePlus;
import ij.gui.Plot;
import ij.measure.ResultsTable;

public class output {

	static void makeResultsTable(ImageDecorrelationAnalysis ida) {
		
		String rtName = "Results";
		
		String units = ida.im.getCalibration().getUnits();
		double pps = ida.im.getCalibration().pixelWidth;
		
		ResultsTable rt;
		if (ij.WindowManager.getWindow(rtName) == null)
			rt = new ResultsTable();
		else 
			rt = ResultsTable.getResultsTable();
		
		// Restrict title to 30 char
		String imTitle = ida.im.getTitle();
		if (imTitle.length() > 30)
			imTitle = imTitle.substring(0, 30);
		
		ImagePlus im = (ImagePlus) ij.WindowManager.getCurrentImage();
		int s = 1; int f = 1; int c = 1;
		if (im != null) {
			s = im.getSlice();
			f = im.getFrame();
			c = im.getChannel();
		}
		
		rt.incrementCounter();
		rt.addLabel(imTitle);
		rt.addValue("Res.", 2*pps/ida.kcMax);
		rt.addValue("units","[" + units + "]");
		rt.addValue("A0", ida.A0);
		rt.addValue("Kc",ida.kcMax);
		rt.addValue("Kc GM", ida.kcGM);
		rt.addValue("Slice", s);
		rt.addValue("Frame", f);
		rt.addValue("Chanel",c);
		
		rt.addValue("rMin", ida.rmin);
		rt.addValue("rMax", ida.rmax);
		rt.addValue("Nr", ida.Nr);
		rt.addValue("Ng",ida.Ng);
		rt.show(rtName);
		
	}
	
	
	static void saveResultsTable(String path) {
		ResultsTable rt = ResultsTable.getResultsTable();
//		IJ.log(path + getID(5) + "_batch_decorrelation_analysis.csv");
		rt.save(path + getID(5) + "_batch_decorrelation_analysis.csv");
	}
	
	static String getID(int n) {
		Calendar cal = Calendar.getInstance();
		String s;
		String y = Integer.toString(cal.get(Calendar.YEAR));
		String m = String.format("%02d", Integer.parseInt(Integer.toString(cal.get(Calendar.MONTH)+1))); // Month starts from 0 to 11
		String d = String.format("%02d", Integer.parseInt(Integer.toString(cal.get(Calendar.DATE))));
		String h = String.format("%02d", Integer.parseInt(Integer.toString(cal.get(Calendar.HOUR))));
		String min = String.format("%02d", Integer.parseInt(Integer.toString(cal.get(Calendar.MINUTE))));
		
		switch (n) {
			case 1 :
				s = y;
				break;
			case 2 :
				s = y + m;
				break;
			case 3 :
				s = y + m + d;
				break;
			case 4 : 
				s = y + m + d + "_" + h;
			case 5 : 
				s = y + m + d + "_" + h + min;
				break;
			default :
				s = y + m + d + "_" + h + min;
				break;
		}
		return s;
	}
	
	static void plotResults(ImageDecorrelationAnalysis ida) {
		
			
			double[] x = new double[ida.Nr];
			
			// Restric title to 30 char
			String imTitle = ida.im.getTitle();
			if (imTitle.length() > 30)
				imTitle = imTitle.substring(0, 30); 
			
			Plot p = new Plot("Decorrelation curve of "+ imTitle, "Normalized spatial frequency", "Decorrelation");

				
			p.setLineWidth(1);
			p.setLimits(0, 1, 0, 1);
			
	// Plot D0
			for (int k = 0; k < ida.d0.length ; k ++)
				x[k] = 0 + (1-0)*(double)k/(double)(ida.Nr-1);

			p.setColor(new Color(0, 0, 0));
//			p.addPoints(x, ida.d0, PlotWindow.LINE);
			p.addPoints(x, ida.d0, Plot.LINE);
			
	// Plot D0 max as a red cross		
			p.setColor(new Color(255, 0, 0));
			double[] ptsX = new double[2];
			double[] ptsY = new double[2];
			ptsX[0] = ida.kc0; ptsX[1] = ida.kc0;
			ptsY[0] = ida.A0; ptsY[1] = ida.A0;
			p.addPoints(ptsX, ptsY, Plot.TRIANGLE);
			
	// Plot high-pass curves 0:Ng-1	
			for (int k = 0; k < ida.d0.length ; k ++)
				x[k] = ida.rmin + (ida.rmax-ida.rmin)*(double)k/(double)(ida.Nr-1);

			double[] dg = new double[ida.Nr];
			for (int k = 0; k < ida.Ng; k++) {
				for (int j = 0; j < ida.Nr; j++)
					dg[j] = ida.d[j][k];
				
				p.setColor(new Color((int)(255*((double)k/((double)ida.Ng-1))), 
							0, (int)(255*(1-(double)k/((double)ida.Ng-1)))));
				p.addPoints(x, dg, Plot.LINE);
			}

	// Plot high-pass curves Ng:2*Ng-1
			for (int k = 0; k < ida.d0.length ; k ++)
				x[k] = ida.rmin2 + (ida.rmax2-ida.rmin2)*(double)k/(double)(ida.Nr-1);

			for (int k = ida.Ng; k < 2*ida.Ng; k++) {
				for (int j = 0; j < ida.Nr; j++)
					dg[j] = ida.d[j][k];
				
				p.setColor(new Color(0, (int)(255*((double)(k-ida.Ng)/((double)ida.Ng-1))),
						(int)(255*((double)(k-ida.Ng)/((double)ida.Ng-1)))));
				p.addPoints(x, dg, Plot.LINE);
			}
			
	// Plot all max
			p.setColor(new Color(0, 255, 0));
			p.addPoints(ida.kc, ida.Ag, Plot.TRIANGLE);
			
	// Plot cutoff frequency as a vertical line
			double[] x0 = new double[2]; double[] y0 = new double[2];
			x0[0] = ida.kcMax; x0[1] = ida.kcMax; 
			y0[0] = 0; y0[1] = 1;
			p.setColor(new Color(0 , 0 ,0));
			p.addPoints(x0, y0, Plot.LINE);

	// Add results on the plot
			double pps = ida.im.getCalibration().pixelWidth;
			String units = ida.im.getCalibration().getUnits();
			
			if (pps == 0)
				pps = 1;
			
			if (ida.kcMax < 0.7) {
				p.addLabel(ida.kcMax + 0.01, 0.1, " Resolution : "+ String.format("%.3f", 2*pps/ida.kcMax) + " ["+units+"]");
				p.addLabel(ida.kcMax + 0.01, 0.2, "kc: "+ String.format("%.3f", ida.kcMax));
				p.addLabel(ida.kcMax + 0.01, 0.3, "A0 : "  + String.format("%.3f", ida.A0));
			}
			else {
				p.addLabel(0.01, 0.1, "      Resolution : "+ String.format("%.3f", 2*pps/ida.kcMax) + " ["+units+"]");
				p.addLabel(0.01, 0.2,  "kc: "+ String.format("%.3f", ida.kcMax));
				p.addLabel(0.01, 0.3,  "A0 : "+ String.format("%.3f", ida.A0));
			}
			
	// Show the plot
			p.show();
		
	};
}
