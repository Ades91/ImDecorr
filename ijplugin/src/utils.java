/* 
% 
%   Set of utilities image processing functions for the ImageDecorrelationAnalysis plugin
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

import ij.ImageStack;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class utils {
	
	
	// This class contains generic utility Image Processing functions required for the analysis

	
	static public FloatProcessor getPreprocessedImage(FloatProcessor im) {
		
		// make it square via zero-pad (mandatory for fht.transform)
		int newSize = Math.max(im.getWidth(), im.getHeight());

		newSize = (int)Math.pow(2, Math.ceil(Math.log((double)newSize)/Math.log(2.0)));
		
		float[] pixelsIn = (float[]) im.getPixels();
		
		// Since the image is already apodized, the padding value is given by any edge pixel
		float meanIm = pixelsIn[0]; 
		
		float[] pixelsOut = new float[newSize*newSize];
		for (int k = 0 ; k < newSize*newSize; k++) { pixelsOut[k] = meanIm; }
		
		int ox = (int) (newSize-im.getWidth())/2;
		int oy = (int) (newSize-im.getHeight())/2;

		int xIn = 0, xOut = ox;
		int yIn = 0, yOut = oy;
		for (xIn = 0, xOut = ox ; xIn <im.getWidth(); xIn++,xOut++)
			for (yIn = 0,yOut = oy; yIn < im.getHeight(); yIn++,yOut++) 
				pixelsOut[yOut*newSize +xOut] = pixelsIn[yIn*im.getWidth() + xIn];
			
		return new FloatProcessor(newSize,newSize,pixelsOut,null);
	}
	
	static public FloatProcessor apodizeEdges(FloatProcessor im) {

		float[] pin = new float[im.getWidth()*im.getHeight()];
		pin = (float[]) im.getPixels();
		
		double dist = 0;
		int off = 20; // number of pixels used for smooth apodization
		double r = im.getWidth()/2-off;
		
		float edgeMean = 0f;
		int count = 0;
		for (int x = 0; x < im.getWidth(); x++)
			for (int y=0;y < im.getHeight(); y++) {
				dist = (x-im.getWidth()/2)*(x-im.getWidth()/2) + (y-im.getHeight()/2)*(y-im.getHeight()/2);
				if (dist > r*r) {
					edgeMean += pin[y*im.getWidth() + x];
					count++;
				}
			}
		edgeMean = edgeMean/count;
		double x0 = 0;
		double y0 = 0;

		for (int x = 0 ; x < im.getWidth(); x++)
			for (int y = 0 ; y < im.getHeight() ; y++) {
			
				x0 = Math.abs(x-im.getWidth()/2);
				y0 = Math.abs(y-im.getHeight()/2);
				if (Math.abs(x0-im.getWidth()/2) <= off | Math.abs(y0-im.getHeight()/2) <= off) { 
					double d = Math.min(Math.abs(x0-im.getWidth()/2), Math.abs(y0-im.getHeight()/2));
					double c = (Math.cos(d*Math.PI/off-Math.PI)+1)/2; 
					
					pin[y*im.getWidth() + x] = (float) (c*(pin[y*im.getWidth() + x]-edgeMean))+edgeMean;
				}
				else if (Math.abs(x-im.getWidth()/2) > im.getWidth()/2 & Math.abs(y-im.getHeight()/2) > im.getHeight()/2) {
					pin[y*im.getWidth() + x] = (float) edgeMean;
				}
			}
		
		FloatProcessor out = new FloatProcessor(im.getWidth(),im.getHeight(),pin);
		
		return out;
	}
	
	static public FloatProcessor substract(FloatProcessor I1, FloatProcessor I2) {
		float[] pixelsI1 = (float[]) I1.getPixels();
		float[] pixelsI2 = (float[]) I2.getPixels();

		float[] pixelsOut = new float[pixelsI1.length];

		for (int k = 0; k < pixelsI1.length; k++)
			pixelsOut[k] = pixelsI1[k]-pixelsI2[k];
		
		return new  FloatProcessor(I1.getWidth(),I1.getHeight(),pixelsOut);
	}
	
	static public FloatProcessor[] normalizeFFT(FloatProcessor[] Ir) {
		
		float[] pixelsInR = (float[]) Ir[0].getPixels();
		float[] pixelsInI = (float[]) Ir[1].getPixels();
		
		FloatProcessor[] out = new FloatProcessor[2];
		out[0] = new FloatProcessor(Ir[0].getWidth(),Ir[0].getHeight());
		out[1] = new FloatProcessor(Ir[0].getWidth(),Ir[0].getHeight());
				
		float[] pixelsOutR = new float[pixelsInR.length];
		float[] pixelsOutI = new float[pixelsInR.length];
		
		for (int k = 0 ; k < pixelsInR.length; k++) {
			float mag = (float) Math.sqrt(pixelsInR[k]*pixelsInR[k] + pixelsInI[k]*pixelsInI[k]);
			if (mag != 0) { // avoid dividing by zero
				pixelsOutR[k] = pixelsInR[k]/mag;
				pixelsOutI[k] = pixelsInI[k]/mag;
			}
		}
		
		out[0].setPixels(pixelsOutR);
		out[1].setPixels(pixelsOutI);
		return out;
	}
	
	static public FloatProcessor[] getComplexFFT(ImageProcessor ip) {
		FHT fht = new FHT(ip);
		fht.transform();

		FloatProcessor[] ret = new FloatProcessor[2];

		ImageStack stack1 = fht.getComplexTransform();
		ret[0] = ((FloatProcessor) stack1.getProcessor(1));
		ret[1] = ((FloatProcessor) stack1.getProcessor(2));

		return ret;
	}
	
}
