/* 
% 
%   ImageDecorrelationAnalysis GUI and input manager
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

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FilenameFilter;

import javax.swing.*;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;


public class ImageDecorrelationAnalysis_ implements PlugIn {

	final private JPanel mainPanel = new JPanel();
	GroupLayout layout = new GroupLayout(mainPanel);
	
	private JPanel settingsPanel;
	private JPanel actionPanel;
	
	private JTextField rminField;
	private JLabel rminLabel;
	private JTextField rmaxField;
	private JLabel rmaxLabel;
	private JTextField NrField;
	private JLabel NrLabel;
	private JTextField NgField;
	private JLabel NgLabel;
	private JCheckBox doPlotCB;
	private JCheckBox batchCB;
	private JCheckBox batchStackCB;
	
	private JButton compute;
	private JButton about;

	
	@Override
	public void run(String arg0) {

// Panels settings
		final JFrame frame;
		frame = new JFrame("Image Decorrelation Analysis");
		mainPanel.setLayout(new BoxLayout(mainPanel,BoxLayout.PAGE_AXIS));
		
		settingsPanel = new JPanel();
		settingsPanel.setLayout(new BoxLayout(settingsPanel,BoxLayout.PAGE_AXIS));
		settingsPanel.setBorder(BorderFactory.createTitledBorder("Settings"));
		settingsPanel.setPreferredSize(new Dimension(300,130));
		
		SpringLayout layout = new SpringLayout();
		actionPanel = new JPanel();
		actionPanel.setLayout(layout);
		actionPanel.setPreferredSize(new Dimension(300,50));

		JPanel subPanel1 = new JPanel();
		subPanel1.setLayout(new FlowLayout());
		JPanel subPanel2 = new JPanel();
		subPanel2.setLayout(new FlowLayout());
		JPanel subPanel3 = new JPanel();
		subPanel2.setLayout(new FlowLayout());
		
// Settings labels and fields 
		rminLabel = new JLabel("Radius min : ");
		rminField = new JTextField("0");
		rminField.setPreferredSize(new Dimension(50,20));
		rminField.setToolTipText("Minimum radius [0,rMax] (normalized frequencies) used for decorrelation analysis.");
		
		rmaxLabel = new JLabel("Radius max : ");
		rmaxField = new JTextField("1");
		rmaxField.setPreferredSize(new Dimension(50,20));
		rmaxField.setToolTipText("Maximum radius [rMin,1] (normalized frequencies) used for decorrelation analysis.");	
		
		NrLabel = new JLabel("Nr : ");
		NrField = new JTextField("50");
		NrField.setPreferredSize(new Dimension(50,20));
		NrField.setToolTipText("[10,100], Sampling of decorrelation curve.");
		
		NgLabel = new JLabel("Ng : ");
		NgField = new JTextField("10");
		NgField.setPreferredSize(new Dimension(50,20));
		NgField.setToolTipText("[5,40], Number of high-pass image analyzed.");
		
		doPlotCB = new JCheckBox();
		doPlotCB.setText("Do plot");
		doPlotCB.setSelected(true);
		doPlotCB.setToolTipText("Plot decorrelation analysis");
		
		batchStackCB = new JCheckBox();
		batchStackCB.setText("Batch stack");
		batchStackCB.setSelected(false);
		batchStackCB.setToolTipText("Batch process all dimensions of the active image");
		
		batchCB = new JCheckBox();
		batchCB.setText("Batch folder");
		batchCB.setSelected(false);
		batchCB.setToolTipText("Batch process all images in folder");
		

		
// Action button settings
		compute = new JButton("Compute");
		compute.setToolTipText("Lauch the decorrelation analysis on the current image.");
		
		compute.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {			
				// grab all parameters and create an ImageDecorrelationAnalysis object
				double rmin,rmax = 0.0; 
				int Nr = 50; int Ng = 10;
				try {
					rmin = Double.parseDouble(rminField.getText());
					rmax = Double.parseDouble(rmaxField.getText());
					Nr = (int)Double.parseDouble(NrField.getText());
					Ng = (int)Double.parseDouble(NgField.getText());
					boolean doPlot = doPlotCB.isSelected();
					boolean batch = batchCB.isSelected();
					boolean batchStack = batchStackCB.isSelected();
					
					// input consistency check
					if (Nr < 10) 	{	Nr = 10; 	NrField.setText(Integer.toString(10));}
					if (Ng < 5) 	{	Ng = 5; 	NgField.setText(Integer.toString(5));}
					if (rmin < 0) 	{	rmin = 0;	rminField.setText(Integer.toString(0));}
					if (rmax > 1) 	{	rmax = 1;	rmaxField.setText(Integer.toString(1));}
					if (rmin == rmax) {
						rmin = 0; rmax = 1;
						rminField.setText(Integer.toString(0));
						rmaxField.setText(Integer.toString(1));
					}
					if (rmin > rmax) {
						rmin = 0; rmax = 1;
						rminField.setText(Integer.toString(0));
						rmaxField.setText(Integer.toString(1));
					}
					
					initDecorrelationAnalysis(rmin,rmax,Nr,Ng,doPlot,batch,batchStack);

						
				} catch (Exception e) {
					IJ.error("Error while initializing the analysis.\n"
							+"Please check input settings.\n" 
							+ e.toString());
				}
			}
		}
		);

		about = new JButton("About");
		about.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) { displayAbout();}
			});
		
// Add all components to their specific panel and set the layout
		actionPanel.add(compute);
		actionPanel.add(about);
		layout.putConstraint(SpringLayout.WEST, compute, 25, SpringLayout.WEST, actionPanel);
		layout.putConstraint(SpringLayout.NORTH, compute, 15, SpringLayout.NORTH, actionPanel);
		layout.putConstraint(SpringLayout.EAST, about, -25, SpringLayout.EAST, actionPanel);
		layout.putConstraint(SpringLayout.NORTH, about, 15, SpringLayout.NORTH, actionPanel);
		subPanel1.add(rminLabel);
		subPanel1.add(rminField);
		subPanel1.add(rmaxLabel);
		subPanel1.add(rmaxField);
		
		subPanel2.add(NrLabel); subPanel2.add(NrField);
		subPanel2.add(NgLabel); subPanel2.add(NgField);
		subPanel3.add(doPlotCB); subPanel3.add(batchStackCB);subPanel3.add(batchCB);
		
		settingsPanel.add(subPanel1);
		settingsPanel.add(subPanel2);
		settingsPanel.add(subPanel3);
		
		mainPanel.add(settingsPanel);
		mainPanel.add(actionPanel);
		
		mainPanel.setLayout(layout);

		frame.add(settingsPanel,BorderLayout.NORTH);
		frame.add(actionPanel,BorderLayout.SOUTH);
		frame.pack();
		frame.setSize(300, 200);
		frame.setLocation(500, 300);
		frame.setResizable(false);
		frame.setVisible(true);
	}
	

	private static void initDecorrelationAnalysis(double rmin,double rmax, int Nr, int Ng, boolean doPlot, boolean batch, boolean batchStack) {
		
		
// File selection management 
		if (batch == false) {
			
			ImagePlus im = (ImagePlus) ij.WindowManager.getCurrentImage();
			
			if (im == null){
				// open file selections tool for images
//				IJ.log("No image found");
				IJ.open();
				im = ij.WindowManager.getCurrentImage();
				im.show();
			}

			// check if current image has a rectangle ROI selection
			if (im.getRoi() != null) {
				if (im.getRoi().getType() == 0) {// Rectangle ROI
//					IJ.log("Number of images : " + im.getStackIndex(im.getNChannels(), im.getNSlices(), im.getNFrames()));
//					IJ.log("IN selection  : # slices " + im.getNSlices() + 
//							", # frames " + im.getNFrames() + ", # channels " + im.getNChannels());
//					im.show();

//					for (int k = 0 ; k < im.getStackIndex(im.getNChannels(), im.getNSlices(), im.getNFrames()); k++) {
//					}
				}
				else {
					IJ.log(   "Non rectangular ROI not supported.\n"
							+ "Deleting ROI and continue analysis.");
					im.deleteRoi();
				}
			}
			
			
			if (batchStack == false){
				// Extract current image
				ImageStack stack = im.getStack();
				ImagePlus imp = new ImagePlus(stack.getSliceLabel(im.getCurrentSlice()),stack.getProcessor(im.getCurrentSlice()));
				imp.setCalibration(im.getCalibration());
				imp.setRoi(im.getRoi());
				imp.setTitle(im.getTitle());
				if(imp.getBitDepth()==24) { // if it is a RGB image
					IJ.run(imp, "8-bit", ""); // Convert to gray scale
				}
				im = imp;
			}
			
			// Image ready to be processed
			
// Create ImageDecorrelationAnalysis object 
			ImageDecorrelationAnalysis IDA = new ImageDecorrelationAnalysis(im,rmin,rmax,Nr,Ng,doPlot,null);
// Run the analysis 
			IDA.startAnalysis();
			
		}
		else { // open selection tool for folder with images in them
			ij.io.DirectoryChooser gd = new ij.io.DirectoryChooser("Please select a directory containing images to be processed.");
			String dirPath = gd.getDirectory();
//			IJ.log("DirPath : " + dirPath);
			File dir = new File(dirPath);
			String[] files = dir.list(new FilenameFilter() {
				public boolean accept(File dir, String name) {
					int index = name.indexOf('.', name.length()-5);
					if (index != -1) {
						String ext = name.substring(index, name.length());
//						IJ.log("files extensions : " + ext);
						if (ext.equals(".tif") || ext.equals(".tiff") || ext.equals(".png") || 
								ext.equals(".ome") || ext.equals(".bmp"))
							return true;
						else 
							return false;
					}
					else 
						return false;
				}
			});

			
			for (int k = 0; k < files.length ; k++) {
				String file = files[k];
//				IJ.log("File # " +  Integer.toString(k) + ": " + file);
				ImagePlus im = IJ.openImage(dirPath + File.separator + file);
				im.show();
				String savePath;
				if (im != null) {
					if (k == files.length-1)
						savePath = dirPath;
					else
						savePath = null;
					
					ImageDecorrelationAnalysis IDA = new ImageDecorrelationAnalysis(im,rmin,rmax,Nr,Ng,doPlot,savePath);
					// Run the analysis
					IDA.startAnalysis();
				}
			}

		}
	}
	
	static void displayAbout() {
		
		/* 
		 * Version 1.1.1: 	Fixed compatibility issues with micro-manager 1.4
		 * Version 1.1.2: 	Removed all logs output. Reordered results table output
		 * Version 1.1.3: 	Added Input variable check
		 * Version 1.1.4: 	Added comment and formating for release
		 * Version 1.1.5: 	Changed resolution criterion to max freq. instead of GM
		 * 					Renamed variables to be consistent with Matlab code and manuscript
		 */
		
		JOptionPane.showMessageDialog(null, 
				  "                             Author : Adrien Descloux\n"
				+ "                            LBEN STI EPFL Switzerland\n"  
				+ "                                           Version 1.1.5\n"
				+"This program is free software: you can redistribute it and/or modify\r\n" + 
				"it under the terms of the GNU General Public License as published by\r\n" 	+ 
				"the Free Software Foundation, either version 2 of the License, or\r\n" 	+ 
				"(at your option) any later version.\r\n" 									+ 
				"This program is distributed in the hope that it will be useful,\r\n" 		+ 
				"but WITHOUT ANY WARRANTY; without even the implied warranty of\r\n" 		+ 
				"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\r\n" 		+ 
				"GNU General Public License for more details.\r\n" 							+ 
				"You should have received a copy of the GNU General Public License\r\n" 	+ 
				"along with this program.  If not, see <http://www.gnu.org/licenses/>", "About this plugin", JOptionPane.INFORMATION_MESSAGE);
	}
	
	public static void main(String[] args) {
		Class<?> clazz = ImageDecorrelationAnalysis_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);
		
		new ij.ImageJ();
		IJ.runPlugIn(clazz.getName(),"");
	}
	
	
}
