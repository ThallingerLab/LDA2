/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2023 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */

package at.tugraz.genome.lda.target.calibration;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.apache.commons.math3.util.Precision;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.entity.ChartEntity;
import org.jfree.chart.entity.XYItemEntity;

import javafx.util.Pair;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class RecalibrationPlotMouseListener implements ChartMouseListener
{
	RecalibrationPlot plot_;
	
	protected RecalibrationPlotMouseListener(RecalibrationPlot plot)
	{
		this.plot_ = plot;
	}

	@Override
	public void chartMouseClicked(ChartMouseEvent e)
	{
		ChartEntity entity = e.getEntity();
		if (entity instanceof XYItemEntity)
		{
			XYItemEntity ce = (XYItemEntity) entity;
			if (!ce.getDataset().equals(plot_.getCurveFitDataset()))
			{
        Double xValue = ce.getDataset().getX(ce.getSeriesIndex(),ce.getItem()).doubleValue();
        Double yValue = ce.getDataset().getY(ce.getSeriesIndex(),ce.getItem()).doubleValue();
        Pair<Double,Double> dataPoint = new Pair<Double,Double>(xValue,yValue);
        
        AnchorPointVO originalAnchorPoint = plot_.getPanel().findLipidClassForDataPoint(dataPoint);
        System.out.println(originalAnchorPoint.getLipidSpecies());
        
        String dataType = CalibrationGraphPanel.PLOT_ALL;
        if (plot_.getDataStandards().contains(originalAnchorPoint))
        {
        	dataType = CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX;
        }
        
        new DataPointDialog(new JFrame(), "Data point selected", originalAnchorPoint, dataType, plot_);
			}
		}
	}

	@Override
	public void chartMouseMoved(ChartMouseEvent e) {}
	
	
	private class DataPointDialog extends JDialog {
		
		private static final long serialVersionUID = 1L;
		
		private RecalibrationPlot plot_;
		private AnchorPointVO dataPoint_;
		

		public DataPointDialog(JFrame parent, String title, AnchorPointVO dataPoint, String dataType, RecalibrationPlot plot) 
		{
	    super(parent, title, true);
	    this.plot_ = plot;
	    this.dataPoint_ = dataPoint;
	    
	    this.setLocationRelativeTo(plot);
	    
	    getContentPane().setLayout(new GridBagLayout());
	    
	    getContentPane().add(new JLabel("Lipid class: "), generateDefaultConstraints(0, 0));
	    getContentPane().add(new JLabel(dataPoint.getLipidClass()), generateDefaultConstraints(7, 0));
	    getContentPane().add(new JLabel("Lipid species: "), generateDefaultConstraints(0, 1));
	    getContentPane().add(new JLabel(dataPoint.getLipidSpecies()), generateDefaultConstraints(7, 1));
	    getContentPane().add(new JLabel("Retention time in original RT-DB: "), generateDefaultConstraints(0, 2));
	    getContentPane().add(new JLabel(String.format("%s min", Precision.round(dataPoint.getxValue(), 2))), generateDefaultConstraints(7, 2));
	    getContentPane().add(new JLabel("Retention time with new experimental conditions: "), generateDefaultConstraints(0, 3));
	    getContentPane().add(new JLabel(String.format("%s min", Precision.round(dataPoint.getyValue()+dataPoint.getxValue(), 2))), generateDefaultConstraints(7, 3));
	    getContentPane().add(new JLabel("Retention time difference: "), generateDefaultConstraints(0, 4));
	    getContentPane().add(new JLabel(String.format("%s min", Precision.round(dataPoint.getyValue(), 2))), generateDefaultConstraints(7, 4));
	    
	    JPanel buttonPane = new JPanel();
	    JButton cancelButton = new JButton("Cancel"); 
	    cancelButton.addActionListener(new ActionListener() {
			  public void actionPerformed(ActionEvent e) 
			  {
			  	cancelButton_actionPerformed(e);
			  }
		  });
	    buttonPane.add(cancelButton, BorderLayout.WEST); 
	    JButton deleteButton = new JButton("Delete this data point");
	    deleteButton.addActionListener(new ActionListener() {
			  public void actionPerformed(ActionEvent e) 
			  {
			  	deleteButton_actionPerformed(e);
			  }
		  });
	    if (dataType.equals(CalibrationFileChooserPanel.DATA_TYPE_STANDARD_MIX))
	    {
	    	deleteButton.setEnabled(false); //standards should not be deletable
	    	deleteButton.setToolTipText("Data points derived from standards cannot be deleted.");
	    }
	    buttonPane.add(deleteButton, BorderLayout.EAST); 
	    getContentPane().add(buttonPane, new GridBagConstraints(0, 5, 5, 1, 0.0, 0.0
	        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 10, 10, 5), 0, 0));
	    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
	    pack(); 
	    setVisible(true);
	  }  
	  
	  private GridBagConstraints generateDefaultConstraints(int x, int y)
	  {
	  	return new GridBagConstraints(x, y, 5, 1, 0.0, 0.0
	        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 10, 10, 5), 0, 0);
	  }
	  
	  private void cancelButton_actionPerformed(ActionEvent e)
	  {
	  	setVisible(false); 
	    dispose();
	  }
	  
	  private void deleteButton_actionPerformed(ActionEvent e)
	  {
	  	this.plot_.getPanel().removeDataPoint(dataPoint_);
	  	setVisible(false); 
	    dispose();
	  }

	}

}
