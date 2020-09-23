/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 16.01.2015
 */
package com.albertozanon.TimeHomogeneousTest;

import static org.junit.Assert.fail;


import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.time.LocalDate;
import java.time.Month;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.Vector;


import org.junit.Assert;
import org.junit.Test;

import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.marketdata.calibration.ParameterObject;
import net.finmath.marketdata.calibration.Solver;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.AnalyticModelFromCurvesAndVols;
import net.finmath.marketdata.model.curves.Curve;
import net.finmath.marketdata.model.curves.CurveInterpolation.ExtrapolationMethod;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationEntity;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationMethod;
import net.finmath.marketdata.model.volatilities.AbstractVolatilitySurface;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterpolation;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveFromDiscountCurve;
import net.finmath.marketdata.products.AnalyticProduct;
import net.finmath.marketdata.products.Swap;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.CalibrationProduct;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromTermStructureModel;
import net.finmath.montecarlo.interestrate.TermStructureModel;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelTimeHomogenousPiecewiseConstant;
import net.finmath.montecarlo.interestrate.models.covariance.TermStructCovarianceModelFromLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.TermStructureCovarianceModelParametric;

import com.albertozanon.TimeHomogeneouLMM.CapletOnCompoundedBackward;
import com.albertozanon.TimeHomogeneouLMM.LIBORMarketModelWithTenorRefinement;

import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.Caplet;
import net.finmath.montecarlo.interestrate.products.Caplet.ValueUnit;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.optimizer.OptimizerFactory;
import net.finmath.optimizer.OptimizerFactoryLevenbergMarquardt;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.Schedule;
import net.finmath.time.ScheduleGenerator;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import net.finmath.time.TimeDiscretizationFromArray.ShortPeriodLocation;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;
import net.finmath.time.daycount.DayCountConvention_ACT_365;

/**
 * This class tests the LIBOR market model and products.
 *
 * @author Christian Fries
 */
public class CalibrationBackwardCapletOnTimeHomogeneous {

	private final int numberOfPaths		= 2000;
	private final int numberOfFactors	= 1;

	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.0000%;-##0.0000%", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterParam		= new DecimalFormat(" #0.00000; -#0.00000", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));

	

	public CalibrationProduct createCalibrationItem( double weight, double maturity, final double moneyness, final double targetVolatility, final ForwardCurve forwardCurve, final DiscountCurve discountCurve) throws CalculationException {
		double strike = 0.004783;
		double dtLibor= 0.5;
		double maturityMinusLengthLibor = maturity - dtLibor;
		CapletOnCompoundedBackward capletBackward = new CapletOnCompoundedBackward(maturityMinusLengthLibor, dtLibor, strike, dtLibor, false);			
		return new CalibrationProduct(capletBackward, targetVolatility, weight);
	}
	
	public static void main(final String[] args) throws CalculationException, SolverException {
		final CalibrationBackwardCapletOnTimeHomogeneous test = new CalibrationBackwardCapletOnTimeHomogeneous();
		test.testATMSwaptionCalibration();
	}

	public void testATMSwaptionCalibration() throws CalculationException, SolverException {
		
		System.out.println("Calibration to Swaptions:");

		final AnalyticModel curveModel = getCalibratedCurve();

		final ForwardCurve forwardCurve = curveModel.getForwardCurve("ForwardCurveFromDiscountCurve(discountCurve-EUR,1D)");

		final DiscountCurve discountCurve = curveModel.getDiscountCurve("discountCurve-EUR");

		final ArrayList<String>			calibrationItemNames	= new ArrayList<>();
		final ArrayList<CalibrationProduct>	calibrationProducts		= new ArrayList<>();

		final String[] atmExpiries = {"1Y", "18M", "2Y", "3Y", "4Y", "5Y", "7Y", "10Y", "15Y", "20Y", "25Y", "30Y" };

		final double[] atmNormalVolatilities = {0.002324,0.002465,0.002790,0.003492,0.004163,0.004634,0.005203,0.00555,0.00547,0.00533};
			
		final LocalDate referenceDate = LocalDate.of(2020, Month.JULY, 31);
		final BusinessdayCalendarExcludingTARGETHolidays cal = new BusinessdayCalendarExcludingTARGETHolidays();
		final DayCountConvention_ACT_365 modelDC = new DayCountConvention_ACT_365();
		for(int i=0; i<atmNormalVolatilities.length; i++ ) 	{
			final LocalDate exerciseDate = cal.getDateFromDateAndOffsetCode(referenceDate, atmExpiries[i]);
			double	exercise		= modelDC.getDaycountFraction(referenceDate, exerciseDate); 
			
			exercise = Math.round(exercise/0.25)*0.25;

			if(exercise < 0.25) {
				continue;
			}
			if(exercise < 1.0) {
				continue;
			}

			final double	moneyness			= 0.0;
			final double	targetVolatility	= atmNormalVolatilities[i];


			final double	weight = 1.0;

			calibrationProducts.add(createCalibrationItem(weight, exercise, moneyness, targetVolatility, forwardCurve, discountCurve));
			calibrationItemNames.add(atmExpiries[i]);
		}
		LIBORModelMonteCarloSimulationModel simulationCalibrated = null;
		final double lastTime	= 21.0;
		final double dt		= 0.005;
		final double dtSimulation	= 0.005;
		System.out.println("dtSimulation: " + dtSimulation);
		final TimeDiscretizationFromArray timeDiscretizationFromArray = new TimeDiscretizationFromArray(0.0, (int) (lastTime / dtSimulation), dtSimulation);
		final TimeDiscretization liborPeriodDiscretization =  new TimeDiscretizationFromArray(0.0, (int) (lastTime / dt), dt);
		int seed=34343;
		System.out.println("Seed: " + seed);
		final BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(timeDiscretizationFromArray, numberOfFactors, numberOfPaths, seed /* seed */);

		//LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTimeHomogenousPiecewiseConstant(timeDiscretizationFromArray, liborPeriodDiscretization, timeToMaturityDiscretization, arrayValues);
	   LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(timeDiscretizationFromArray, liborPeriodDiscretization, 0.2/100.0, 0.1/100.0, 0.15, 0.05/100.0, true); //0.20/100.0, 0.05/100.0, 0.10, 0.05/100.0, 

		final LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, 0.05, false);

		//AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelExponentialForm5Param(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, new double[] { 0.20/100.0, 0.1/100.0, 0.10, 0.05/100.0, 0.10} );
		AbstractLIBORCovarianceModelParametric	covarianceModelParametric = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretizationFromArray, liborPeriodDiscretization, volatilityModel, correlationModel);
		TermStructureCovarianceModelParametric termStructureCovarianceModel = new TermStructCovarianceModelFromLIBORCovarianceModelParametric(null, covarianceModelParametric );

		
		int k=0;
		// Set model properties
		final Map<String, Object> properties = new HashMap<>();

		final Double accuracy = 1E-12;
		final int maxIterations = 400;
		final int numberOfThreads = 6;
		final OptimizerFactory optimizerFactory = new OptimizerFactoryLevenbergMarquardt(maxIterations, accuracy, numberOfThreads);

		final double[] parameterStandardDeviation = new double[termStructureCovarianceModel.getParameter().length];
		final double[] parameterLowerBound = new double[termStructureCovarianceModel.getParameter().length];
		final double[] parameterUpperBound = new double[termStructureCovarianceModel.getParameter().length];
		Arrays.fill(parameterStandardDeviation, k==0 ? 0.0020/100.0 : 0.2/100.0);
		Arrays.fill(parameterLowerBound, Double.NEGATIVE_INFINITY);
		Arrays.fill(parameterUpperBound, Double.POSITIVE_INFINITY);

		//				optimizerFactory = new OptimizerFactoryCMAES(accuracy, maxIterations, parameterLowerBound, parameterUpperBound, parameterStandardDeviation);

		// Set calibration properties (should use our brownianMotion for calibration - needed to have to right correlation).
		final Map<String, Object> calibrationParameters = new HashMap<>();
		calibrationParameters.put("accuracy", accuracy);
		calibrationParameters.put("brownianMotion", brownianMotion);
		calibrationParameters.put("parameterStep", k == 0 ? new Double(1E-6) : new Double(5E-5) );
		calibrationParameters.put("optimizerFactory", optimizerFactory);
		properties.put("calibrationParameters", calibrationParameters);
		
		
		final CalibrationProduct[] calibrationItemsLMM = new CalibrationProduct[calibrationItemNames.size()];
		for(int i=0; i<calibrationItemNames.size(); i++) {
			calibrationItemsLMM[i] = new CalibrationProduct(calibrationProducts.get(i).getProduct(),calibrationProducts.get(i).getTargetValue(),calibrationProducts.get(i).getWeight());
		}
		
		final TimeDiscretization liborPeriodDiscretizationDaily = new TimeDiscretizationFromArray(0.0, 21.0, 0.005, ShortPeriodLocation.SHORT_PERIOD_AT_START);
		final TimeDiscretization liborPeriodDiscretizationWeekly = new TimeDiscretizationFromArray(0.0, 21.0, 0.025, ShortPeriodLocation.SHORT_PERIOD_AT_START);
		final TimeDiscretization liborPeriodDiscretizationMonthly = new TimeDiscretizationFromArray(0.0, 21.0, 0.125, ShortPeriodLocation.SHORT_PERIOD_AT_START);
		final TimeDiscretization liborPeriodDiscretizationQuarterly = new TimeDiscretizationFromArray(0.0, 21.0, 0.25, ShortPeriodLocation.SHORT_PERIOD_AT_START);
		final TimeDiscretization liborPeriodDiscretizationSemiannual = new TimeDiscretizationFromArray(0.0, 21.0, 0.5, ShortPeriodLocation.SHORT_PERIOD_AT_START);

		final TermStructureModel liborMarketModelCalibrated = new LIBORMarketModelWithTenorRefinement(
				new TimeDiscretization[] { liborPeriodDiscretizationDaily, liborPeriodDiscretizationWeekly, liborPeriodDiscretizationMonthly,liborPeriodDiscretizationQuarterly,liborPeriodDiscretizationSemiannual},
				new Integer[] { 5, 4, 3, 2, 200 },
				curveModel,
				forwardCurve,
				new DiscountCurveFromForwardCurve(forwardCurve),
				termStructureCovarianceModel,
				calibrationItemsLMM,
				properties		
		);


		System.out.println("\nCalibrated parameters are:");
		final double[] param = ((LIBORMarketModelWithTenorRefinement) liborMarketModelCalibrated).getCovarianceModel().getParameter();
		//		((AbstractLIBORCovarianceModelParametric) liborMarketModelCalibrated.getCovarianceModel()).setParameter(param);
		for (final double p : param) {
			System.out.println(p);
		}

		final EulerSchemeFromProcessModel process = new EulerSchemeFromProcessModel(liborMarketModelCalibrated,brownianMotion);
		simulationCalibrated = new LIBORMonteCarloSimulationFromTermStructureModel( process);
	
		System.out.println("\nValuation on calibrated model:");
		double deviationSum			= 0.0;
		double deviationSquaredSum	= 0.0;
		
		for (int i = 0; i < calibrationProducts.size(); i++) {
			final AbstractLIBORMonteCarloProduct calibrationProduct = calibrationProducts.get(i).getProduct();
			try {
				final double valueModel = calibrationProduct.getValue(simulationCalibrated);
				final double valueTarget = calibrationProducts.get(i).getTargetValue().getAverage();
				final double error = valueModel-valueTarget;
				deviationSum += error;
				deviationSquaredSum += error*error;
				System.out.println(calibrationItemNames.get(i) + "\t" + "Model: " + formatterValue.format(valueModel) + "\t Target: " + formatterValue.format(valueTarget) + "\t Deviation: " + formatterDeviation.format(valueModel-valueTarget));// + "\t" + calibrationProduct.toString());
			}
			catch(final Exception e) {
			}
		}

		final double averageDeviation = deviationSum/calibrationProducts.size();
		System.out.println("Mean Deviation:" + formatterValue.format(averageDeviation));
		System.out.println("RMS Error.....:" + formatterValue.format(Math.sqrt(deviationSquaredSum/calibrationProducts.size())));
		System.out.println("__________________________________________________________________________________________\n");

		Assert.assertTrue(Math.abs(averageDeviation) < 1E-2);
	

//		// CAPLET ON BACKWARD LOOKING RATE SEMESTRALI
		DecimalFormat formatterTimeValue = new DecimalFormat("##0.00;");
		DecimalFormat formatterVolValue = new DecimalFormat("##0.00000;");
		DecimalFormat formatterAnalytic = new DecimalFormat("##0.000;");
		DecimalFormat formatterPercentage = new DecimalFormat(" ##0.000%;-##0.000%", new DecimalFormatSymbols(Locale.ENGLISH));
	
		double strike = 0.004783;
		double dtLibor = 0.5;
		
		double[] mktData = new double[] {/* 6M 0.00167, */ /* 12M*/ 0.00201, /* 18M*/ 0.00228, /* 2Y */ 0.00264, 0.0, /* 3Y */ 0.0033, /* 4Y */0.00406, /* 5Y */ 0.00455, /* 6Y - NA */ 0.0, /* 7Y */0.00513, /* 8Y- NA */0.0, /* 9Y */0.0, /* 10Y */0.00550,0.0,0.0,0.0,0.0, /* 15Y */0.00544,0.0,0.0,0.0,0.0, /* 20Y */0.0053,0.0,0.0,0.0};
	
		int mktDataIndex = 0;

		
		//Results with CALIBRATED model
		System.out.println("\n results on CALIBRATED model \n");
		double beginLiborTime = 0.5;
		while(beginLiborTime < liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1)) {
			Caplet capletCassical = new Caplet(beginLiborTime, dtLibor, strike, dtLibor, false, ValueUnit.NORMALVOLATILITY);
			CapletOnCompoundedBackward backwardCaplet = new CapletOnCompoundedBackward(beginLiborTime, dtLibor, strike, dtLibor,false);
			double capletMaturity = beginLiborTime+dtLibor;
			double impliedVolClassical = capletCassical.getValue(simulationCalibrated);
			double impliedVolBackward = backwardCaplet.getValue(simulationCalibrated);
			double analyticFormulaPaper = Math.sqrt(1+0.5/(beginLiborTime*3));
			double ratioImpliedVol = impliedVolBackward/impliedVolClassical;
			double error = (analyticFormulaPaper-ratioImpliedVol)/analyticFormulaPaper;
			
			if (beginLiborTime<2.5) {		//da i valori del caplet per maturity 1.5Y, 2Y, 2.5Y,.	
				if (mktData[mktDataIndex] == 0.0) {
				System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
						+ formatterTimeValue.format(capletMaturity) + "; "
						+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
						+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
						+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
						+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
						+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
					+ formatterPercentage.format(error) );
				beginLiborTime+=0.5;
				mktDataIndex+=1;
				}
				else {

					double ratioMktVol =impliedVolBackward/mktData[mktDataIndex];
					double errorMkt = (analyticFormulaPaper-ratioMktVol)/analyticFormulaPaper;

					System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
							+ formatterTimeValue.format(capletMaturity) + "; "
							+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
							+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
							+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
							+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
							+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
							+ formatterPercentage.format(error) 
							+  "\t" +"Market Caplet Vol. " + formatterPercentage.format(mktData[mktDataIndex])
							+  "\t" +"Market Ratio: "	+ formatterAnalytic.format(ratioMktVol)	
							+  "\t" +"Market Error: "	+ formatterPercentage.format(errorMkt)
							);
					beginLiborTime+=0.5;
					mktDataIndex+=1;
				}
				
			}
			else {	//secondo loop da i valori del caplet per maturity 4Y,5Y,...,21
	
				if (mktData[mktDataIndex] == 0.0) {
					System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
							+ formatterTimeValue.format(capletMaturity) + "; "
							+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
							+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
							+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
							+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
							+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
						+ formatterPercentage.format(error) );
					beginLiborTime+=1;
					mktDataIndex+=1;
					}
					else {
						double ratioMktVol =impliedVolBackward/mktData[mktDataIndex];
						double errorMkt = (analyticFormulaPaper-ratioMktVol)/analyticFormulaPaper;

						System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
								+ formatterTimeValue.format(capletMaturity) + "; "
								+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
								+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
								+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
								+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
								+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
								+ formatterPercentage.format(error) 
								+  "\t" +"Market Caplet Vol. " + formatterPercentage.format(mktData[mktDataIndex])
								+  "\t" +"Market Ratio: "	+ formatterAnalytic.format(ratioMktVol)	
								+  "\t" +"Market Error: "	+ formatterPercentage.format(errorMkt)
								);
						beginLiborTime+=1;
						mktDataIndex+=1;
					}
			}
		}
		
		System.out.println("\n Backward looking rate:");
		for(beginLiborTime=0.0; beginLiborTime<liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1); beginLiborTime+=1) {
			double endLiborTime=beginLiborTime+dtLibor;
			RandomVariable backwardLookingRate =  simulationCalibrated.getBackward(endLiborTime, beginLiborTime,endLiborTime);
			double avgBackwardLookingRate =backwardLookingRate.getAverage();
			System.out.println("Backward B(" + formatterTimeValue.format(beginLiborTime) +  ", " + formatterTimeValue.format(endLiborTime) + "; "  + formatterTimeValue.format(endLiborTime)+ ")."+ "\t" + "Average: " + avgBackwardLookingRate);
		}
		
		System.out.println("\n LIBOR rate:");
		for(beginLiborTime=0.0; beginLiborTime<liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1); beginLiborTime+=1) {
			double endLiborTime=beginLiborTime+dtLibor;
			RandomVariable libor =  simulationCalibrated.getLIBOR(beginLiborTime, beginLiborTime,endLiborTime);
			double avgBackwardLookingRate = libor.getAverage();
			System.out.println("Backward B(" + formatterTimeValue.format(beginLiborTime) +  ", " + formatterTimeValue.format(endLiborTime) + "; "  + formatterTimeValue.format(beginLiborTime)+ ")."+ "\t" + "Average: " + avgBackwardLookingRate);
		}
		
		
		
		final CalibrationProduct[] calibrationItems = new CalibrationProduct[0];

		final TermStructureModel  liborMarketModelNOTCalibrated = new LIBORMarketModelWithTenorRefinement (
				new TimeDiscretization[] { liborPeriodDiscretizationDaily, liborPeriodDiscretizationWeekly, liborPeriodDiscretizationMonthly,liborPeriodDiscretizationQuarterly,liborPeriodDiscretizationSemiannual},
				new Integer[] { 5, 4, 3, 2, 200 },
				//new TimeDiscretization[] { liborPeriodDiscretizationDaily,liborPeriodDiscretizationSemiannual},
				//new Integer[] { 2, 200 },
				curveModel,
				forwardCurve,
				new DiscountCurveFromForwardCurve(forwardCurve),
				termStructureCovarianceModel,
				calibrationItems,
				null		
		);

		final EulerSchemeFromProcessModel process2 = new EulerSchemeFromProcessModel(liborMarketModelNOTCalibrated,brownianMotion);
		LIBORModelMonteCarloSimulationModel  simulationNONCalibrated = new LIBORMonteCarloSimulationFromTermStructureModel ( process2);

		
		System.out.println("\n results on NON CALIBRATED model \n");
		
		 beginLiborTime = 0.5;
		mktDataIndex = 0;

		while(beginLiborTime < liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1)) {
			Caplet capletCassical = new Caplet(beginLiborTime, dtLibor, strike, dtLibor, false, ValueUnit.NORMALVOLATILITY);
			CapletOnCompoundedBackward backwardCaplet = new CapletOnCompoundedBackward(beginLiborTime, dtLibor, strike, dtLibor,false);
			double capletMaturity = beginLiborTime+dtLibor;
			double impliedVolClassical = capletCassical.getValue(simulationNONCalibrated);
			double impliedVolBackward = backwardCaplet.getValue(simulationNONCalibrated);
			double analyticFormulaPaper = Math.sqrt(1+0.5/(beginLiborTime*3));
			double ratioImpliedVol = impliedVolBackward/impliedVolClassical;
			double error = (analyticFormulaPaper-ratioImpliedVol)/analyticFormulaPaper;

			
			if (beginLiborTime<2.5) {		//da i valori del caplet per maturity 1.5Y, 2Y, 2.5Y,.	
				if (mktData[mktDataIndex] == 0.0) {
				System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
						+ formatterTimeValue.format(capletMaturity) + "; "
						+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
						+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
						+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
						+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
						+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
					+ formatterPercentage.format(error) );
				beginLiborTime+=0.5;
				mktDataIndex+=1;
				}
				else {

					double ratioMktVol =impliedVolBackward/mktData[mktDataIndex];
					double errorMkt = (analyticFormulaPaper-ratioMktVol)/analyticFormulaPaper;

					System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
							+ formatterTimeValue.format(capletMaturity) + "; "
							+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
							+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
							+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
							+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
							+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
							+ formatterPercentage.format(error) 
							+  "\t" +"Market Caplet Vol. " + formatterPercentage.format(mktData[mktDataIndex])
							+  "\t" +"Market Ratio: "	+ formatterAnalytic.format(ratioMktVol)	
							+  "\t" +"Market Error: "	+ formatterPercentage.format(errorMkt)
							);
					beginLiborTime+=0.5;
					mktDataIndex+=1;
				}
				
			}
			else {	//secondo loop da i valori del caplet per maturity 4Y,5Y,...,21
	
				if (mktData[mktDataIndex] == 0.0) {
					System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
							+ formatterTimeValue.format(capletMaturity) + "; "
							+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
							+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
							+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
							+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
							+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
						+ formatterPercentage.format(error) );
					beginLiborTime+=1;
					mktDataIndex+=1;
					}
					else {
						double ratioMktVol =impliedVolBackward/mktData[mktDataIndex];
						double errorMkt = (analyticFormulaPaper-ratioMktVol)/analyticFormulaPaper;

						System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
								+ formatterTimeValue.format(capletMaturity) + "; "
								+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
								+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
								+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
								+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
								+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
								+ formatterPercentage.format(error) 
								+  "\t" +"Market Caplet Vol. " + formatterPercentage.format(mktData[mktDataIndex])
								+  "\t" +"Market Ratio: "	+ formatterAnalytic.format(ratioMktVol)	
								+  "\t" +"Market Error: "	+ formatterPercentage.format(errorMkt)
								);
						beginLiborTime+=1;
						mktDataIndex+=1;
					}
			}
		}	
		System.out.println("\n Backward looking rate:");
		for(beginLiborTime=0.0; beginLiborTime<liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1); beginLiborTime+=1) {
			double endLiborTime=beginLiborTime+dtLibor;
			RandomVariable backwardLookingRate =  simulationNONCalibrated.getBackward(endLiborTime, beginLiborTime,endLiborTime);
			double avgBackwardLookingRate =backwardLookingRate.getAverage();
			System.out.println("Backward B(" + formatterTimeValue.format(beginLiborTime) +  ", " + formatterTimeValue.format(endLiborTime) + "; "  + formatterTimeValue.format(endLiborTime)+ ")."+ "\t" + "Average: " + avgBackwardLookingRate);
		}
		
		System.out.println("\n LIBOR rate:");
		for(beginLiborTime=0.0; beginLiborTime<liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1); beginLiborTime+=1) {
			double endLiborTime=beginLiborTime+dtLibor;
			RandomVariable libor =  simulationNONCalibrated.getLIBOR(beginLiborTime, beginLiborTime,endLiborTime);
			double avgBackwardLookingRate = libor.getAverage();
			System.out.println("Backward B(" + formatterTimeValue.format(beginLiborTime) +  ", " + formatterTimeValue.format(endLiborTime) + "; "  + formatterTimeValue.format(beginLiborTime)+ ")."+ "\t" + "Average: " + avgBackwardLookingRate);
		}
	
	
	
	
	}
		
	
	
	
	public AnalyticModel getCalibratedCurve() throws SolverException {
			final String[] maturity					= { "1D", "7D", "14D", "21D", "1M", "2M", "3M", "4M", "5M", "6M", "7M", "8M", "9M", "12M", "15M", "18M", "21M", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "11Y", "12Y", "15Y", "20Y", "25Y", "30Y", "40Y", "50Y" };
			final String[] frequency				= { "tenor", "tenor", "tenor",  "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual"};
			final String[] frequencyFloat			= { "tenor", "tenor", "tenor",  "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual"};
			final String[] daycountConventions	    = { "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360"};
			final String[] daycountConventionsFloat	= { "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360"};
			final double[] rates					= { -0.0055, -0.00553, -0.00553, -0.00553, -0.00553, -0.00555, -0.00556, -0.00559, -0.00564, -0.00568, -0.00572, -0.00577, -0.00581, -0.00592, -0.00601, -0.00608, -0.00613, -0.00619, -0.00627, -0.00622, -0.00606, -0.00582, -0.00553, -0.00519, -0.00482, -0.00442, -0.00402, -0.00362, -0.00261, -0.00189, -0.00197, -0.0023, -0.00286, -0.00333};
			final HashMap<String, Object> parameters = new HashMap<>();

			parameters.put("referenceDate", LocalDate.of(2020, Month.JULY, 31));
			parameters.put("currency", "EUR");
			parameters.put("forwardCurveTenor", "1D");
			parameters.put("maturities", maturity);
			parameters.put("fixLegFrequencies", frequency);
			parameters.put("floatLegFrequencies", frequencyFloat);
			parameters.put("fixLegDaycountConventions", daycountConventions);
			parameters.put("floatLegDaycountConventions", daycountConventionsFloat);
			parameters.put("rates", rates);

		return getCalibratedCurve(null, parameters);
	}

	private static AnalyticModel getCalibratedCurve(final AnalyticModel model2, final Map<String, Object> parameters) throws SolverException {

		final LocalDate	referenceDate		= (LocalDate) parameters.get("referenceDate");
		final String	currency			= (String) parameters.get("currency");
		final String	forwardCurveTenor	= (String) parameters.get("forwardCurveTenor");
		final String[]	maturities			= (String[]) parameters.get("maturities");
		final String[]	frequency			= (String[]) parameters.get("fixLegFrequencies");
		final String[]	frequencyFloat		= (String[]) parameters.get("floatLegFrequencies");
		final String[]	daycountConventions	= (String[]) parameters.get("fixLegDaycountConventions");
		final String[]	daycountConventionsFloat	= (String[]) parameters.get("floatLegDaycountConventions");
		final double[]	rates						= (double[]) parameters.get("rates");

		Assert.assertEquals(maturities.length, frequency.length);
		Assert.assertEquals(maturities.length, daycountConventions.length);
		Assert.assertEquals(maturities.length, rates.length);

		Assert.assertEquals(frequency.length, frequencyFloat.length);
		Assert.assertEquals(daycountConventions.length, daycountConventionsFloat.length);

		final int		spotOffsetDays = 2;
		final String	forwardStartPeriod = "0D";

		final String curveNameDiscount = "discountCurve-" + currency;

		/*
		 * We create a forward curve by referencing the same discount curve, since
		 * this is a single curve setup.
		 *
		 * Note that using an independent NSS forward curve with its own NSS parameters
		 * would result in a problem where both, the forward curve and the discount curve
		 * have free parameters.
		 */
		final ForwardCurve forwardCurve		= new ForwardCurveFromDiscountCurve(curveNameDiscount, referenceDate, forwardCurveTenor);

		// Create a collection of objective functions (calibration products)
		final Vector<AnalyticProduct> calibrationProducts = new Vector<>();
		final double[] curveMaturities	= new double[rates.length+1];
		final double[] curveValue			= new double[rates.length+1];
		final boolean[] curveIsParameter	= new boolean[rates.length+1];
		curveMaturities[0] = 0.0;
		curveValue[0] = 1.0;
		curveIsParameter[0] = false;
		for(int i=0; i<rates.length; i++) {
//			scheduleRec=scheduleReceiveLeg; schedulePay=schedulePayLeg
			final Schedule schedulePay = ScheduleGenerator.createScheduleFromConventions(referenceDate, spotOffsetDays, forwardStartPeriod, maturities[i], frequency[i], daycountConventions[i], "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), -2, 0);
			final Schedule scheduleRec = ScheduleGenerator.createScheduleFromConventions(referenceDate, spotOffsetDays, forwardStartPeriod, maturities[i], frequencyFloat[i], daycountConventionsFloat[i], "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), -2, 0);

			curveMaturities[i+1] = Math.max(schedulePay.getPayment(schedulePay.getNumberOfPeriods()-1),scheduleRec.getPayment(scheduleRec.getNumberOfPeriods()-1));
			curveValue[i+1] = 1.0;
			curveIsParameter[i+1] = true;
			calibrationProducts.add(new Swap(schedulePay, null, rates[i], curveNameDiscount, scheduleRec, forwardCurve.getName(), 0.0, curveNameDiscount));
		} 

		final InterpolationMethod interpolationMethod = InterpolationMethod.LINEAR;

		// Create a discount curve
		final DiscountCurveInterpolation discountCurveInterpolation = DiscountCurveInterpolation.createDiscountCurveFromDiscountFactors(
				curveNameDiscount								/* name */,
				referenceDate	/* referenceDate */,
				curveMaturities	/* maturities */,
				curveValue		/* discount factors */,
				curveIsParameter,
				interpolationMethod ,
				ExtrapolationMethod.CONSTANT,
				InterpolationEntity.LOG_OF_VALUE
				);

		/*
		 * Model consists of the two curves, but only one of them provides free parameters.
		 */
		AnalyticModel model = new AnalyticModelFromCurvesAndVols(new Curve[] { discountCurveInterpolation, forwardCurve });

		/*
		 * Create a collection of curves to calibrate
		 */
		final Set<ParameterObject> curvesToCalibrate = new HashSet<>();
		curvesToCalibrate.add(discountCurveInterpolation);

		/*
		 * Calibrate the curve
		 */
		final Solver solver = new Solver(model, calibrationProducts);
		final AnalyticModel calibratedModel = solver.getCalibratedModel(curvesToCalibrate);
		System.out.println("Solver reported acccurary....: " + solver.getAccuracy());

		Assert.assertEquals("Calibration accurarcy", 0.0, solver.getAccuracy(), 1E-3);

		// Get best parameters
		final double[] parametersBest = calibratedModel.getDiscountCurve(discountCurveInterpolation.getName()).getParameter();

		// Test calibration
		model			= calibratedModel;

//--> praticamente il calibrato valore di calibrationProducts (che sono swap) deve essere = 0 --> penso perchè il prezzo dello swap deve essere = 0 per definizione, quindi in teoria gli interest definiti sono L_i - S e quando calibri questa differenza deve essere =0 (e se così è allora il prezzo dello swap è 0) --> avrebbe anche senso perchè se noti il payer ha tassi prima negativi e poi positivi, mentre quell'altro ha tassi sempre = 0 (cioè L=S, cioè è quello che paga il fixed swap rate)
		double squaredErrorSum = 0.0;
		for(final AnalyticProduct c : calibrationProducts) {
			final double value = c.getValue(0.0, model);
			final double valueTaget = 0.0;
			final double error = value - valueTaget;
			squaredErrorSum += error*error;
		}
		final double rms = Math.sqrt(squaredErrorSum/calibrationProducts.size());

		System.out.println("Independent checked acccurary: " + rms);

		System.out.println("Calibrated discount curve: ");
		for(int i=0; i<curveMaturities.length; i++) {
			final double maturity = curveMaturities[i];
			System.out.println(maturity + "\t" + calibratedModel.getDiscountCurve(discountCurveInterpolation.getName()).getDiscountFactor(maturity));
		}
		return model;
	}

}
