/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 16.01.2015
 */
package com.albertozanon.TimeHomogeneousTest;

import static org.junit.Assert.fail;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
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
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;

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
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.marketdata.products.AnalyticProduct;
import net.finmath.marketdata.products.Swap;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.models.HullWhiteModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModelWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelWithTenorRefinement;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.BlendedLocalVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.DisplacedLocalVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecayWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelExponentialForm5Param;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFourParameterExponentialFormWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelPiecewiseConstant;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelTimeHomogenousPiecewiseConstant;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification;
import net.finmath.montecarlo.interestrate.models.covariance.ShortRateVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.ShortRateVolatilityModelAsGiven;
import net.finmath.montecarlo.interestrate.models.covariance.TermStructCovarianceModelFromLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.TermStructureCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.TermStructureTenorTimeScalingInterface;
import net.finmath.montecarlo.interestrate.models.covariance.TermStructureTenorTimeScalingPicewiseConstant;
import net.finmath.montecarlo.interestrate.models.covariance.VolatilityReductionMercurioModel;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.Caplet;
import net.finmath.montecarlo.interestrate.products.CapletOnBackwardLookingRate;
import net.finmath.montecarlo.interestrate.products.SwaptionSimple;
import net.finmath.montecarlo.interestrate.products.Caplet.ValueUnit;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.Optimizer;
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
public class TestCaplet {

	private final int numberOfPaths		= 10000;
	private final int numberOfFactors	= 1;

	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.0000%;-##0.0000%", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterParam		= new DecimalFormat(" #0.00000; -#0.00000", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));

	

	public CalibrationProduct createCalibrationItem( double weight, double maturityMinusLengthLibor, final double moneyness, final double targetVolatility, final ForwardCurve forwardCurve, final DiscountCurve discountCurve) throws CalculationException {
		double strike = 0.004783;
		double dtLibor= 0.5;
		CapletOnBackwardLookingRate capletBackward = new CapletOnBackwardLookingRate(maturityMinusLengthLibor, dtLibor, strike, dtLibor, false);			
		return new CalibrationProduct(capletBackward, targetVolatility, weight);
	}
		
    @Test
	public void testATMSwaptionCalibration() throws CalculationException, SolverException {
		
		 final RandomVariableFactory randomVariableFactory = new RandomVariableFromArrayFactory();

		System.out.println("Calibration to Swaptions:");

		final AnalyticModel curveModel = getCalibratedCurve();

		final ForwardCurve forwardCurve = curveModel.getForwardCurve("ForwardCurveFromDiscountCurve(discountCurve-EUR,1D)");

		final DiscountCurve discountCurve = curveModel.getDiscountCurve("discountCurve-EUR");

		final ArrayList<String>			calibrationItemNames	= new ArrayList<>();
		final ArrayList<CalibrationProduct>	calibrationProducts		= new ArrayList<>();

		final String[] atmExpiries = {"1Y", "18M", "2Y", "3Y", "4Y", "5Y", "7Y", "10Y", "15Y", "20Y", "25Y", "30Y" };

		final double[] atmNormalVolatilities = {0.002324,0.002465,0.002790,0.003492,0.004163,0.004634,0.005203,0.00555,0.00547,0.00533,0.00513,0.0049};
			
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
		final double lastTime	= 31.0;
		final double dtLibor	= 0.5;
		final double dt	= 0.125;
		final TimeDiscretization timeDiscretizationFromArray = new TimeDiscretizationFromArray(0.0, (int) (lastTime / dt), dt);
		final TimeDiscretization liborPeriodDiscretization = new TimeDiscretizationFromArray(0.0, (int) (lastTime / dtLibor), dtLibor);
		
		final BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(timeDiscretizationFromArray, numberOfFactors, numberOfPaths, 210996 /* seed */);
//		TimeDiscretization optionMaturityDiscretization = new TimeDiscretizationFromArray(0.00, 1.00, 2.00, 3.00, 4.00, 5.00,  7.00,  10.0, 15.0, 20.0, 25.0, 31.0);
		TimeDiscretization timeToMaturityDiscretization = new TimeDiscretizationFromArray(0.00, 1.00, 2.00, 3.00, 4.00, 5.00,  7.00,  10.0, 15.0, 20.0, 25.0, 31.0);		// needed if you use LIBORVolatilityModelPiecewiseConstantWithMercurioModification: TimeDiscretization  = new TimeDiscretizationFromArray(0.0, 0.25, 0.50, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0);
	//	TimeDiscretization timeToMaturityDiscretization = new TimeDiscretizationFromArray(0.00, 0.25, 0.5, 0.75, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9,00, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 31.0);
		//TimeDiscretization timeToMaturityDiscretization = new TimeDiscretizationFromArray(0.00, 31, 1.0);

		double[] arrayValues = new double [timeToMaturityDiscretization.getNumberOfTimes()];
		for (int i=0; i<timeToMaturityDiscretization.getNumberOfTimes(); i++) {arrayValues[i]= 0.2/100;}

		//final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization, timeToMaturityDiscretization, arrayValues);
		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialFormWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization, 0.002, 0.0005, 0.10, 0.0005, true);
		//final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelPiecewiseConstantWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization,optionMaturityDiscretization,timeToMaturityDiscretization, 0.50 / 100);
		
		final LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecayWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, 0.05, false);
		final AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretizationFromArray, liborPeriodDiscretization, volatilityModel, correlationModel);

		// Create blended local volatility model with fixed parameter (0=lognormal, > 1 = almost a normal model).
		final AbstractLIBORCovarianceModelParametric covarianceModelDisplaced = new DisplacedLocalVolatilityModel(covarianceModelParametric, 1.0/0.25, false /* isCalibrateable */);
		final AbstractLIBORCovarianceModelParametric covarianceModelReducedVolatility = new VolatilityReductionMercurioModel(covarianceModelDisplaced);
		
		final Map<String, Object> properties = new HashMap<>();
		System.out.println("Number of volatility parameters: " + volatilityModel.getParameter().length);

		properties.put("measure", LIBORMarketModelFromCovarianceModelWithMercurioModification.Measure.SPOT.name());
		properties.put("stateSpace", LIBORMarketModelFromCovarianceModelWithMercurioModification.StateSpace.NORMAL.name());

		final Double accuracy = new Double(1E-8);	
		final int maxIterations = 200;
		final int numberOfThreads = 1;
		final OptimizerFactory optimizerFactory = new OptimizerFactoryLevenbergMarquardt(maxIterations, accuracy, numberOfThreads);

		final double[] parameterStandardDeviation = new double[covarianceModelParametric.getParameterAsDouble().length];
		final double[] parameterLowerBound = new double[covarianceModelParametric.getParameterAsDouble().length];
		final double[] parameterUpperBound = new double[covarianceModelParametric.getParameterAsDouble().length];
		Arrays.fill(parameterStandardDeviation, 0.20/100.0);
		Arrays.fill(parameterLowerBound, 0.0);
		Arrays.fill(parameterUpperBound, Double.POSITIVE_INFINITY);

		final Map<String, Object> calibrationParameters = new HashMap<>();
		calibrationParameters.put("accuracy", accuracy);
		calibrationParameters.put("brownianMotion", brownianMotion);
		calibrationParameters.put("optimizerFactory", optimizerFactory);
		calibrationParameters.put("parameterStep", new Double(1E-4));
		properties.put("calibrationParameters", calibrationParameters);

		final long millisCalibrationStart = System.currentTimeMillis();

		final CalibrationProduct[] calibrationItemsLMM = new CalibrationProduct[calibrationItemNames.size()];
		for(int i=0; i<calibrationItemNames.size(); i++) {
			calibrationItemsLMM[i] = new CalibrationProduct(calibrationProducts.get(i).getProduct(),calibrationProducts.get(i).getTargetValue(),calibrationProducts.get(i).getWeight());
		}
		final LIBORMarketModel mercurioModelCalibrated = LIBORMarketModelFromCovarianceModelWithMercurioModification.of(
				liborPeriodDiscretization,
				curveModel,
				forwardCurve,
				new DiscountCurveFromForwardCurve(forwardCurve),
				randomVariableFactory,
				covarianceModelReducedVolatility,
				calibrationItemsLMM, properties);

		final long millisCalibrationEnd = System.currentTimeMillis();
		
//-------------------------------------------------------------------------------- fine calibrazione volatility------------------------------
		System.out.println("\nCalibrated parameters are:");
		final double[] param = ((AbstractLIBORCovarianceModelParametric)((LIBORMarketModelFromCovarianceModelWithMercurioModification) mercurioModelCalibrated).getCovarianceModel()).getParameterAsDouble();
		for (final double p : param) {
			System.out.println(p);
		}

		final EulerSchemeFromProcessModel process = new EulerSchemeFromProcessModel(brownianMotion);
		final LIBORModelMonteCarloSimulationModel simulationMercurioCalibrated = new LIBORMonteCarloSimulationFromLIBORModel(mercurioModelCalibrated, process);

		System.out.println("\nValuation on calibrated model:");
		double deviationSum			= 0.0;
		double deviationSquaredSum	= 0.0;
		for (int i = 0; i < calibrationProducts.size(); i++) {
			final AbstractLIBORMonteCarloProduct calibrationProduct = calibrationProducts.get(i).getProduct();
			try {
				final double valueModel = calibrationProduct.getValue(simulationMercurioCalibrated);
				final double valueTarget = calibrationProducts.get(i).getTargetValue().getAverage();
				final double error = valueModel-valueTarget;
				deviationSum += error;
				deviationSquaredSum += error*error;
				System.out.println(calibrationItemNames.get(i) + "\t" + "Model: " + formatterValue.format(valueModel) + "\t Target: " + formatterValue.format(valueTarget) + "\t Deviation: " + formatterDeviation.format(valueModel-valueTarget));// + "\t" + calibrationProduct.toString());
			}
				catch(final Exception e) {
			}
		}
		System.out.println("Time required for calibration of volatilities...: " + (millisCalibrationEnd-millisCalibrationStart)/1000.0 + " s.");

		final double averageDeviation = deviationSum/calibrationProducts.size();
		System.out.println("Mean Deviation:" + formatterValue.format(averageDeviation));
		System.out.println("RMS Error.....:" + formatterValue.format(Math.sqrt(deviationSquaredSum/calibrationProducts.size())));
		System.out.println("__________________________________________________________________________________________\n");
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
