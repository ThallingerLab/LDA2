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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.util.ArrayList;

import javax.swing.JPanel;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.function.Function2D;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import org.jfree.chart.util.ShapeUtils;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class RecalibrationPlot extends JPanel
{ 
	private static final long serialVersionUID = 1L;
	
	ArrayList<AnchorPointVO> dataAll_;
	ArrayList<AnchorPointVO> dataStandards_;
	CalibrationGraphPanel panel_;
	XYPlot plot_;
	XYDataset curveFitDataset_;
	XYDataset scatterPlotStandardsDataset_;
	XYDataset scatterPlotDataset_;


	public RecalibrationPlot(ArrayList<AnchorPointVO> data, ArrayList<AnchorPointVO> dataStandards, 
			RecalibrationRegression regression, Dimension dimension, CalibrationGraphPanel panel, XYPlot previousPlot)
  {
		this.dataAll_ = new ArrayList<AnchorPointVO>(data);
		this.dataStandards_ = dataStandards;
		this.panel_ = panel;
		this.dataAll_.removeAll(this.dataStandards_);
		
		plot_ = new XYPlot();
    Font fontTitle = new Font("Dialog", Font.BOLD, 25);
    Font fontLabel = new Font("Dialog", Font.PLAIN, 20);
    Font tickLabel = new Font("Dialog", Font.PLAIN, 14);
    
    plot_.setDomainAxis(0, getAxis("Retention time original RT-DB /min", fontLabel, tickLabel));
    plot_.setRangeAxis(0, getAxis("Retention time difference original vs new RT-DB /min", fontLabel, tickLabel));
    
    if (previousPlot != null)
    {
    	copyAxisRanges(plot_, previousPlot);
    }
    
    plot_.mapDatasetToDomainAxis(0, 0);
    plot_.mapDatasetToRangeAxis(0, 0);
    
    this.curveFitDataset_ = getCurveFitDataset(regression);
    plot_.setDataset(0, this.curveFitDataset_);
    plot_.setRenderer(0, getCurveRenderer(fontLabel));
    
    this.scatterPlotStandardsDataset_ = getScatterPlotDataset(this.dataStandards_, String.format("Standards (%s calibrants)", this.dataStandards_.size()));
    plot_.setDataset(1, this.scatterPlotStandardsDataset_);
    plot_.setRenderer(1, getScatterRenderer(fontLabel, new Color(204, 0, 0)));
    
    this.scatterPlotDataset_ = getScatterPlotDataset(this.dataAll_, String.format("Other (%s calibrants)", this.dataAll_.size()));
    plot_.setDataset(2, this.scatterPlotDataset_);
    plot_.setRenderer(2, getScatterRenderer(fontLabel, new Color(153, 204, 255)));

    JFreeChart chart = new JFreeChart("RT-DB mapping",
    		fontTitle, plot_, true);
    ChartPanel chartPanel = new ChartPanel(chart);
    chartPanel.addChartMouseListener(new RecalibrationPlotMouseListener(this));
    chartPanel.setPreferredSize(dimension);
    
    this.add(chartPanel);
  }
	
	private void copyAxisRanges(XYPlot plot, XYPlot previousPlot)
	{
		ValueAxis domainAxis = plot.getDomainAxis();
		ValueAxis previousDomainAxis = previousPlot.getDomainAxis();
		ValueAxis rangeAxis = plot.getRangeAxis();
		ValueAxis previousRangeAxis = previousPlot.getRangeAxis();
		
		if (!previousDomainAxis.isAutoRange())
		{
			domainAxis.setAutoRange(false);
			domainAxis.setRange(previousDomainAxis.getRange());
		}
		
		if (!previousRangeAxis.isAutoRange())
		{
			rangeAxis.setAutoRange(false);
			rangeAxis.setRange(previousRangeAxis.getRange());
		}
	}
	
	private NumberAxis getAxis(String label, Font fontLabel, Font tickLabel)
	{
		NumberAxis axis = new NumberAxis(label);
	  axis.setLabelFont(fontLabel);
	  axis.setTickLabelFont(tickLabel);
	  return axis;
	}
	
	private XYLineAndShapeRenderer getCurveRenderer(Font fontLabel)
	{
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
    renderer.setSeriesPaint(0, new Color(0, 0, 102));
    renderer.setSeriesPaint(1, new Color(0, 0, 102));
    renderer.setDefaultLegendTextFont(fontLabel);
    renderer.setAutoPopulateSeriesStroke(false);
    renderer.setDefaultStroke(new BasicStroke(2f));
    return renderer;
	}
  
  private XYDataset getCurveFitDataset(RecalibrationRegression regression)
  {
  	int num = 10000;
  	XYSeriesCollection dataset = new XYSeriesCollection();
  	XYSeries dataSeries = new XYSeries("Fit");
  	double upperLimit = regression.getUpperRTLimit();
  	double lowerLimit = regression.getLowerRTLimit();
  	double step = (upperLimit-lowerLimit)/num;
  	for (int i=0; i<=num;i++)
  	{
  		double x = lowerLimit+i*step;
  		if (regression.getFunction().isValidPoint(x))
  		{
  			double y = regression.getFunction().value(x);
    		dataSeries.add(x,y);
  		}
  	}
  	dataset.addSeries(dataSeries);
  	return dataset;
  }

  public static XYSeries sampleFunctionOverY(Function2D f, double start,
      double end, int samples, Comparable<?> seriesKey)
  {
    XYSeries series = new XYSeries(seriesKey, false);
    double step = (end - start) / (samples - 1);
    for (int i = 0; i < samples; i++)
    {
      double y = start + step * i;
      series.add(f.getValue(y), y);
    }
    return series;
  }
  
  private XYLineAndShapeRenderer getScatterRenderer(Font fontLabel, Color fillPaint)
	{
  	XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(false, true);
    renderer.setDefaultLegendTextFont(fontLabel);
    renderer.setDefaultFillPaint(fillPaint);
    renderer.setDefaultOutlinePaint(new Color(64, 64, 64));
    renderer.setUseFillPaint(true);
    renderer.setUseOutlinePaint(true);
    renderer.setSeriesShape(0, ShapeUtils.createDiamond(5));
    return renderer;
	}
  
  private XYDataset getScatterPlotDataset(ArrayList<AnchorPointVO> data, String label)
  {
    XYSeriesCollection dataset = new XYSeriesCollection();
    XYSeries dataSeries = new XYSeries(label);
    for (AnchorPointVO dataPoint : data)
    {
    	dataSeries.add(dataPoint.getxValue(), dataPoint.getyValue());
    }
    dataset.addSeries(dataSeries);
    return dataset;
  }

	public ArrayList<AnchorPointVO> getDataAll()
	{
		return dataAll_;
	}

	public ArrayList<AnchorPointVO> getDataStandards()
	{
		return dataStandards_;
	}

	public CalibrationGraphPanel getPanel()
	{
		return panel_;
	}

	public XYDataset getCurveFitDataset()
	{
		return curveFitDataset_;
	}

	public XYDataset getScatterPlotStandardsDataset()
	{
		return scatterPlotStandardsDataset_;
	}

	public XYDataset getScatterPlotDataset()
	{
		return scatterPlotDataset_;
	}
	
	public XYPlot getXYPlot()
	{
		return plot_;
	}
  
  
}
