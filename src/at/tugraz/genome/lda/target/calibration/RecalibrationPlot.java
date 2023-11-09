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

import at.tugraz.genome.lda.utils.Pair;


public class RecalibrationPlot extends JPanel
{ 
	private static final long serialVersionUID = 1L;
	
	ArrayList<Pair<Double,Double>> dataAll_;
	ArrayList<Pair<Double,Double>> dataStandards_;
	CalibrationGraphPanel panel_;
	XYPlot plot_;
	XYDataset curveFitDataset_;
	XYDataset scatterPlotStandardsDataset_;
	XYDataset scatterPlotDataset_;


	public RecalibrationPlot(ArrayList<Pair<Double,Double>> data, ArrayList<Pair<Double,Double>> dataStandards, 
			RecalibrationRegression regression, Dimension dimension, CalibrationGraphPanel panel, XYPlot previousPlot)
  {
		this.dataAll_ = new ArrayList<Pair<Double,Double>>(data);
		this.dataStandards_ = dataStandards;
		this.panel_ = panel;
		this.dataAll_.removeAll(this.dataStandards_);
		
		plot_ = new XYPlot();
    Font fontTitle = new Font("Dialog", Font.BOLD, 25);
    Font fontLabel = new Font("Dialog", Font.PLAIN, 20);
    Font tickLabel = new Font("Dialog", Font.PLAIN, 14);
    
    plot_.setDomainAxis(0, getAxis("Retention time original target list /min", fontLabel, tickLabel));
    plot_.setRangeAxis(0, getAxis("Retention time difference original vs new target list /min", fontLabel, tickLabel));
    
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

    JFreeChart chart = new JFreeChart("Target list recalibration curve",
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
  	XYSeries dataSeries = new XYSeries("Fitted recalibration curve");
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
  
  private XYDataset getScatterPlotDataset(ArrayList<Pair<Double,Double>> data, String label)
  {
    XYSeriesCollection dataset = new XYSeriesCollection();
    XYSeries dataSeries = new XYSeries(label);
    for (Pair<Double,Double> dataPoint : data)
    {
    	dataSeries.add(dataPoint.getKey(), dataPoint.getValue());
    }
    dataset.addSeries(dataSeries);
    return dataset;
  }

	public ArrayList<Pair<Double,Double>> getDataAll()
	{
		return dataAll_;
	}

	public ArrayList<Pair<Double,Double>> getDataStandards()
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
