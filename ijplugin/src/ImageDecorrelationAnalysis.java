/* 
% 
%   Main ImageDecorrelationAnalysis class
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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class ImageDecorrelationAnalysis{

	public ImagePlus im;
	public double rmin;
	public double rmax;
	public double rmin2;
	public double rmax2;
	
	public int Nr;
	public int Ng;
	public boolean doPlot;
	public double[] d0;
	public double[][] d;
	public double kcGM,kcMax,kc0 = 0;
	public double AGM,AMax,A0 = 0;
	public double[] kc;
	public double[] Ag;
	public ImageProcessor imRef;
	public ImageProcessor imN;
	public int s,f,c;
	public String savePath;
	
	ImageDecorrelationAnalysis (ImagePlus im, double rmin, double rmax, int Nr, int Ng, boolean doPlot, String savePath) {
		this.im = im;
		this.rmin = rmin; 	this.rmax = rmax;
		this.rmin2 = rmin; 	this.rmax2 = rmax;
		this.Nr = Nr;	 	this.Ng = Ng;
		this.doPlot = doPlot;
		this.d0 = new double[this.Nr];
		this.d = new double[this.Nr][2*this.Ng];
		this.kc = new double[2*this.Ng];
		this.Ag = new double[2*this.Ng];
		this.kcGM = 0; this.AGM = 0;
		this.kcMax = 0; this.AMax = 0;
		this.s = 1; this.f = 1; this.c = 1;
		this.savePath = savePath;
	}
	

	public void startAnalysis() {
		Thread t1 = new Thread(new Runnable() {
			public void run() {
				runAnalysis();
			}
			
		});
		t1.start();
		try {
			t1.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		saveResults();
	}
	
	public void runAnalysis() {
		
//		this.im = this.im.crop();
		int nS = this.im.getNSlices();
		int nF = this.im.getNFrames();
		int nC = this.im.getNChannels();
//		IJ.log("Processing image : # slices " + nS + ", # frames " + nF + ", # channels " + nC );
		
		// loop through all images 
		for (int s = 1; s <= nS; s++)
			for (int f = 1 ; f <= nF; f++)
				for (int c = 1; c <= nC ; c++){
					
//					IJ.log("Processing image : slice " + s + ", frame " + f + ", channel " + c);
					this.c = c; this.s = s; this.f = f;
					this.im.setPosition(c, s, f);
					ImageStack stack = this.im.getStack();

					ImagePlus imp = new ImagePlus(stack.getSliceLabel(this.im.getCurrentSlice()),stack.getProcessor(this.im.getCurrentSlice()));
					imp.setCalibration(this.im.getCalibration());
					imp.setRoi(this.im.getRoi());
					imp.setTitle(this.im.getTitle());

//					IJ.showProgress(0);
					
// Grab current image and crop if possible
					this.imRef = imp.getProcessor().crop();
					
//					IJ.log("Input image preprocessing (W : " + imRef.getWidth() + ", H : " + imRef.getHeight() + ")");
					
// pre-processing steps : Apodization, make square image with a 2^N size
					// Make FloatProcessor
					this.imRef = this.imRef.toFloat(0, null);
					FloatProcessor imf = new FloatProcessor(this.imRef.getWidth(),this.imRef.getHeight(),(float[]) this.imRef.getPixels());
					imf = utils.apodizeEdges(imf);
					
// make image square with a 2^N size
					FloatProcessor temp = utils.getPreprocessedImage(imf);
					this.imRef = temp;
//					IJ.log("Compute Fourier transform inverse");
//					FloatProcessor[] im2 = getComplexiFFT(temp);
//					ImagePlus I = new ImagePlus();
//					I.setProcessor(this.imRef);
//					I.show();
					
//					IJ.log("GetComplexFFT");
// Compute Fourier transfrom of pre-processed image
					FloatProcessor[] Ir = utils.getComplexFFT(temp);
					
					// Remove mean from Ir
					Ir[0].set(Ir[0].getWidth()/2, Ir[0].getHeight()/2, 0);
					Ir[1].set(Ir[0].getWidth()/2, Ir[0].getHeight()/2, 0);
					
//					ImagePlus I = new ImagePlus();
//					I.setProcessor(Ir[0]);
//					I.show();
					
// Compute d0 curve
//					long t0 = System.currentTimeMillis();
//					IJ.log("Compute Dcor curve");
					computeD0(Ir); 
//					long t1 = System.currentTimeMillis();
//					IJ.log("Computation time for a single curve : " + (t1-t0));
					
// Compute d0 max position and amplitude
					double[] out = getDcorrMax(this.d0,0,1);
					this.kc0 = out[0];
					this.A0 = out[1];
					
//					IJ.log("D0 end : " + this.d0[this.Nr-1]);
//					IJ.log("A0 : "     + this.A0);
//					IJ.log("res0 : "   + this.res0);
					
// Refine the analysis with high-pass filtering
//					t0 = System.currentTimeMillis();
					computeD(this.imRef); 
//					IJ.log("");
//					t1 = System.currentTimeMillis();
//					IJ.log("Computation time for all the curves : " + (t1-t0));		
									
// Plot results
					if (this.doPlot) {
//						IJ.log("Plot results");
						output.plotResults(this);
					}
					
// Display results in table
					output.makeResultsTable(this);
					
				}
		
		saveResults();
		
	}
	
	private void saveResults() {
		if (this.savePath != null)
			output.saveResultsTable(this.savePath);
	}
	
	
	private void computeD(ImageProcessor Ip) {
		
		double gMax;
		int count = 0;
		
		if (this.kc0 == 0)
			gMax = Ip.getWidth()/2;
		else
			gMax = 2/this.kc0;

		double gMin = 0.14;
		GaussianBlur gb = new GaussianBlur();
		
		FloatProcessor imRef = Ip.convertToFloatProcessor();
		FloatProcessor imBlur = Ip.convertToFloatProcessor();
		
		imRef = utils.getPreprocessedImage(imRef);
		
		FloatProcessor[] Ir = utils.getComplexFFT(imRef);
		FloatProcessor[] I = utils.normalizeFFT(Ir);
		
		double crmin = this.rmin;
		double crmax = this.rmax;
		boolean[] mask0 = getMask(Ir[0].getWidth(),1);
//		IJ.log("start loop");
		
		for (int refine = 0; refine < 2; refine++) {
			for (int k = 0; k < this.Ng; k++) {
				// High-pass filtering of Ip
				imBlur = imRef.convertToFloatProcessor();

				double sig =  Math.exp(Math.log(gMin) + (Math.log(gMax)-Math.log(gMin))*((double)k/((double)this.Ng-1)));
				gb.blurGaussian(imBlur,sig,sig,0.01);

				imBlur = utils.substract(imRef,imBlur);

				Ir = utils.getComplexFFT(imBlur);
				I = utils.normalizeFFT(Ir);

//				boolean[] mask = new boolean[Ir[0].getWidth()*Ir[0].getHeight()];
				float[] pixelsIrR = (float[]) Ir[0].getPixels();
				float[] pixelsIrI = (float[]) Ir[1].getPixels();
				
				float[] pixelsIR = (float[]) I[0].getPixels();
				float[] pixelsII = (float[]) I[1].getPixels();
				double cr = getCorrCoefNorm(pixelsIrR,pixelsIrI,mask0);
				
				for (int i = 0; i <this.Nr; i++) {
					double r = crmin + (crmax-crmin)*(double)i/(double)(this.Nr-1);
//					mask = getMask(Ir[0].getWidth(),r*r);

					this.d[i][count] = getCorrCoef(pixelsIrR,pixelsIrI,pixelsIR,pixelsII,r,cr);
//					IJ.showProgress(0.1 + 0.9*(((double)count2)/((this.Nr-1)*2*(this.Ng-1))));

//					count2++;
				}
				count++;
				
			}
			
			if (refine == 0) {
				// first set of curves are computed => refine high-pass and sampling
				double[] kc = new double[this.Ng+1];
				double[] A  = new double[this.Ng+1];
				double[] dg = new double[this.Nr];
				double[] result = new double[2];
				
				for (int j = 0 ; j < this.Ng ; j++) {
					for (int h = 0 ; h < this.Nr ; h++)
						dg[h] = this.d[h][j];
					
					result = getDcorrMax(dg,crmin,crmax);
//					IJ.log("dg curve kc : " + result[0] + ", A : " + result[1] );
					kc[j] = result[0];
					A[j] = result[1];
					this.kc[j] = result[0];
					this.Ag[j] = result[1];
				}
				
				// Add d0 as potential candidate for best fit
				result = getDcorrMax(this.d0,crmin,crmax);
//				IJ.log("D0 curve kc : " + result[0] + ", A : " + result[1] );
				kc[this.Ng] = result[0];
				A[this.Ng] = result[1];
				this.kc[this.Ng] = result[0];
				this.Ag[this.Ng] = result[1];
				
				// process peaks positions with geometric mean
				double[] resultsGM = getBestScore(kc,A);
				double[] resultsMax = getMaxScore(kc,A);
				
				// refine the analysis according to results
				crmin = Math.min(resultsGM[0],resultsMax[0])-0.05;
				if (crmin < this.rmin)
					crmin = this.rmin;
				
				crmax = Math.max(resultsGM[0],resultsMax[0])+0.3;
				if(crmax > this.rmax)
					crmax = this.rmax;
				
				//	Store refined data for final results and plotting
				this.rmin2 = crmin;
				this.rmax2 = crmax;
				
				double ind1 = Math.min(resultsGM[2],resultsMax[2])-1;
				double ind2 = Math.max(resultsGM[2],resultsMax[2]);
				if (ind2 < this.Ng) {
					double gTemp = Math.exp(Math.log(gMin) + (Math.log(gMax)-Math.log(gMin))*(ind1/((double)this.Ng-1)));
					gMax = Math.exp(Math.log(gMin) + (Math.log(gMax)-Math.log(gMin))*(ind2/((double)this.Ng-1)));
					gMin = gTemp;
				}
				else { // D0 is the best fit
					gMax = gMin;
					gMin = 2/Ip.getWidth();
				}

			}
			else {
				// Compute other half of of curves
				double[] kc = new double[this.Ng];
				double[] A = new double[this.Ng];
				double[] dg = new double[this.Nr];
				double[] result = new double[2];
				
				for (int j = 0 ; j < this.Ng ; j++) {
					for (int h = 0 ; h < this.Nr ; h++)
						dg[h] = this.d[h][j+this.Ng];
					
					result = getDcorrMax(dg,crmin,crmax);
					kc[j] = result[0];
					A[j] = result[1];
					this.kc[j+this.Ng] = result[0];
					this.Ag[j+this.Ng] = result[1];
				}
				
				// process peaks positions with geometric mean
				double[] resultsGM = getBestScore(kc,A);
				this.kcGM = resultsGM[0];
				this.AGM = resultsGM[1];
				// process peaks for maximum frequency
				double[] resultsMax = getMaxScore(kc,A);
				this.kcMax = resultsMax[0];
				this.AMax = resultsMax[1];
			}
//			IJ.log("End of loop");
		}
	}
	
	private double[] getBestScore(double[] kc, double[] A) {

		double[] gm = new double[kc.length];
		double[] out = new double[3];
		double gmMax = 0;
		for (int k = 0 ; k < kc.length ; k++) {
			gm[k] = kc[k]*A[k];
			if(gm[k] > gmMax) {
				gmMax = gm[k];
				out[0] = kc[k];
				out[1] = A[k];
				out[2] = k;
			}	
		}
		return out;
	}
	
	private double[] getMaxScore(double[] kc, double[] A) {
		double out[] = new double[3];
		double kcmax = 0;
		for (int k = 0; k < kc.length; k++) {
			if (kc[k] > kcmax) {
				kcmax = kc[k];
				out[0] = kc[k];
				out[1] = A[k];
				out[2] = k;
			}
		}
		return out;
	}
	
	
	private void computeD0(FloatProcessor[] Ir) {
		
		// Compute normalized version of Ir
//		long t1 = System.currentTimeMillis();
		FloatProcessor[] I = utils.normalizeFFT(Ir);
//		IJ.log("Normalize FFT time computation time : " + (System.currentTimeMillis()-t1));
//		ImagePlus Ip = new ImagePlus();
//		Ip.setProcessor(I[0]);
//		Ip.show();
		
//		boolean[] mask = new boolean[Ir[0].getWidth()*Ir[0].getHeight()];
		boolean mask0[] = getMask(Ir[0].getWidth(),1);
		
		float[] pixelsIrR = (float[]) Ir[0].getPixels();
		float[] pixelsIrI = (float[]) Ir[1].getPixels();
		
		float[] pixelsIR = (float[]) I[0].getPixels();
		float[] pixelsII = (float[]) I[1].getPixels();

		// Precompute the normalization correlation coeficient for Ir
		double cr = getCorrCoefNorm(pixelsIrR,pixelsIrI,mask0);
		
		for (int k = 0; k <this.Nr; k++) {
//			t1 = System.currentTimeMillis();
			// d0 is always computed on the full frequency range
			double r = 0 + (1-0)*(double)k/(double)(this.Nr-1);
//			mask = getMask(Ir[0].getWidth(),r*r);

			this.d0[k] = getCorrCoef(pixelsIrR,pixelsIrI,pixelsIR,pixelsII,r,cr);
//			IJ.log("Computation of a single cc time computation time : " + (System.currentTimeMillis()-t1));
			IJ.showProgress(0.1*(k/(this.Nr-1)));
		}
	}
	
	private double[] getDcorrMax(double[] d, double r1, double r2) {

		double[] t = d.clone();
		double[] out  = getMax(d,0,this.Nr);
		double tempMin[] = getMin(d,0,this.Nr);
		int dLength = t.length;
		double dt = 0.001;
		
		while(out[0] == dLength-1) { // while the maximum of the curve is the last point of the curve
			t[dLength-1] = 0;
			dLength -= 1;
			if (dLength == 0) {
				out[0] = 0;
				out[1] = 0;
				break;
			}
			else {
				out = getMax(t,0,this.Nr);
				tempMin = getMin(t,(int)out[0],dLength-1);
				if (t[(int)out[0]] - tempMin[1] > dt) {
					break;
				}
				else {
					t[(int)out[0]] = tempMin[1];
					out[0] = dLength-1;
				}
			}
		}
				
		out[0] = r1 + (r2-r1)*out[0]/(this.Nr-1); // convert max position into normalized frequency
		return out;
	}

	private double[] getMax(double[] array, int x1, int x2) {
		double[] out = new double[2];
		
		out[0] = x1;
		out[1] = array[x1];
		for (int k = x1;k < x2; k++) {
			if (array[k] > out[1]) {
				out[0] = k;
				out[1] = array[k];
			}	
		}
		return out;
	}
	
	private double[] getMin(double[] array, int x1, int x2) {
		double[] out = new double[2];
		
		out[0] = x1;
		out[1] = array[x1];
		for (int k = x1;k < x2; k++) {
			if (array[k] < out[1]) {
				out[0] = k;
				out[1] = array[k];
			}	
		}
		return out;
	}
	
	private double getCorrCoef(float[] pixelsIrR, float[] pixelsIrI, float[] pixelsIR, float[] pixelsII, double r, double cr) {
		
		double d = 0;
		double c = 0;
		double r2 = r*r;
		double dist = 0;
		int W = this.imRef.getWidth();
		int H = this.imRef.getHeight();
		int k = 0;
		
		if (r != 0){
			int ox =  (int)(W*(1-r)/2);
			int oy =  (int)(H*(1-r)/2);
			int w = (int)(W*r);
			int h = (int)(H*r);
			
			for (int xx = ox; xx < w+ox; xx ++)
				for (int yy = oy; yy < h+oy; yy++) {
					dist = (xx-W/2)*(xx-W/2) + (yy-H/2)*(yy-H/2);
					dist = 4*dist/(W*W);
					if (dist < r2) {
						k = xx*this.imRef.getWidth() + yy;
						d  += pixelsIrR[k]*pixelsIR[k] + pixelsIrI[k]*pixelsII[k];
						c  += pixelsIR[k]*pixelsIR[k] + pixelsII[k]*pixelsII[k];
					}
				}
					
			d = d/(cr*Math.sqrt(c));
		}
		
		return d;
	}
	
	
	boolean[] getMask(int w,double r2) {
		
		r2 = r2*w*w/4;
		boolean[] mask = new boolean[w*w];
		for (int k = 0;k < w*w; k++) {
			int x = (k/w);
			int y = k-w*x;
			x -= w/2;
			y -= w/2;
			
			double dist = x*x  + y*y;
			if (dist < r2)
				mask[k] = true;
			else
				mask[k] = false;
		}
		return mask;
	}
	

	
	private double getCorrCoefNorm(float[] pixelsInR, float[] pixelsInI, boolean[] mask) {
		double c = 0;
		
		for (int k = 0; k < pixelsInR.length; k++) {
			if (mask[k])
				c += pixelsInR[k]*pixelsInR[k] + pixelsInI[k]*pixelsInI[k];
		}
		
		return Math.sqrt(c);
	}


}
