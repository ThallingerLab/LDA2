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

package at.tugraz.genome.lda.target.experiment;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JPanel;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.function.Function2D;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class TotalIsotopeEffectPlot extends JPanel
{ 
	private static final long serialVersionUID = 1L;


	public TotalIsotopeEffectPlot(Vector<MatchedPartnerVO> data, IsotopeEffectRegression regression, Dimension dimension)
  {
    XYPlot plot = new XYPlot();
    Font fontTitle = new Font("Dialog", Font.BOLD, 25);
    Font fontLabel = new Font("Dialog", Font.PLAIN, 20);
    Font tickLabel = new Font("Dialog", Font.PLAIN, 14);
    int domainAxisMax = regression.getMaxNumDeuteriumAllowed();
//    Stroke dashed =  new BasicStroke(0.5f);
//    Stroke dashed =  new BasicStroke(1.0f,BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] {10.0f}, 0.0f);
    XYDataset scatterPlotDataset = getScatterPlotDataset(data);
    plot.setDataset(0, scatterPlotDataset);
    XYLineAndShapeRenderer scatterRenderer = new XYLineAndShapeRenderer(false, true);
    scatterRenderer.setDefaultLegendTextFont(fontLabel);
    plot.setRenderer(0, scatterRenderer);
    NumberAxis domainAxis = new NumberAxis("Number of deuterium atoms");
    domainAxis.setLabelFont(fontLabel);
    domainAxis.setTickLabelFont(tickLabel);
    domainAxis.setRange(0, domainAxisMax);
    plot.setDomainAxis(0, domainAxis);
    NumberAxis rangeAxis = new NumberAxis("Total Isotope Effect (TIE)");
    rangeAxis.setLabelFont(fontLabel);
    rangeAxis.setTickLabelFont(tickLabel);
    rangeAxis.setRange(1, findMaximumForRangeAxis(data, regression, domainAxisMax));
    plot.setRangeAxis(0, rangeAxis);
    plot.mapDatasetToDomainAxis(0, 0);
    plot.mapDatasetToRangeAxis(0, 0);

//    XYDataset functionDataset = 
//        getFunctionDataset(0.8, 0.5, 1.2, minY, maxY);
    XYDataset curveFitDataset = getCurveFitDataset(regression);
    plot.setDataset(1, curveFitDataset);
    XYLineAndShapeRenderer curveRenderer = new XYLineAndShapeRenderer(true, false);
    curveRenderer.setDefaultLegendTextFont(fontLabel);
    curveRenderer.setSeriesFillPaint(1, new Color(0xFF, 0xFF, 0xFF));
//    curveRenderer.setAutoPopulateSeriesStroke(false);
//    curveRenderer.setSeriesStroke(1,new BasicStroke(0.0005f));
    plot.setRenderer(1, curveRenderer);

    JFreeChart chart = new JFreeChart("Isotope effect on retention time",
    		fontTitle, plot, true);
    ChartPanel chartPanel = new ChartPanel(chart);
    chartPanel.setPreferredSize(dimension);
    
    this.add(chartPanel);
  }
  
  private double findMaximumForRangeAxis(Vector<MatchedPartnerVO> data, IsotopeEffectRegression regression, int domainAxisMax)
  {
  	double max = 1.0;
  	for (MatchedPartnerVO matched : data)
    {
    	if (!matched.isUseForCalibration()) continue;
			double standardRT = matched.getStandard().getPreciseRT();
			double isotopologueRT = matched.getIsotopologue().getPreciseRT();
			double totalIsotopeEffect = standardRT / isotopologueRT;
			max = totalIsotopeEffect > max ? totalIsotopeEffect : max;
    }
  	max = regression.getIsotopeEffect(domainAxisMax) > max ? regression.getIsotopeEffect(domainAxisMax) : max;
  	return (max - 1.0)/10.0 + max; //add one tenth to the top
  }
  
  private XYDataset getCurveFitDataset(IsotopeEffectRegression regression)
  {
  	int num = 10000;
  	XYSeriesCollection dataset = new XYSeriesCollection();
  	XYSeries dataSeries = new XYSeries("Fit");
  	double step = new Double(regression.getMaxNumDeuteriumAllowed())/num;
  	for (int i=0; i<=num;i++)
  	{
  		double x = i*step;
  		double y = regression.getIsotopeEffect(x);
  		dataSeries.add(x,y);
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


  private XYDataset getScatterPlotDataset(Vector<MatchedPartnerVO> data)
  {
      XYSeriesCollection dataset = new XYSeriesCollection();
      XYSeries dataSeries = new XYSeries("Matched isotopologues");
      for (MatchedPartnerVO matched : data)
      {
      	if (!matched.isUseForCalibration()) continue;
  			double standardRT = matched.getStandard().getPreciseRT();
  			double isotopologueRT = matched.getIsotopologue().getPreciseRT();
  			try
  			{
  				Hashtable<String,Integer> chemicalFormula = StaticUtils.categorizeFormula(matched.getIsotopologue().getChemicalFormula());
  				int numberDeuterium = chemicalFormula.get("D");
  				double totalIsotopeEffect = standardRT / isotopologueRT;
  				dataSeries.add(numberDeuterium, totalIsotopeEffect);
  			}
  			catch (ChemicalFormulaException ex)
  			{
  				ex.printStackTrace();
  			}
      }
      dataset.addSeries(dataSeries);
      return dataset;
  }
}
