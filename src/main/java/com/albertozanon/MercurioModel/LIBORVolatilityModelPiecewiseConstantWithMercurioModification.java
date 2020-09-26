/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 08.08.2005
 */
package com.albertozanon.MercurioModel;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;

public class LIBORVolatilityModelPiecewiseConstantWithMercurioModification extends LIBORVolatilityModel {

	private static final long serialVersionUID = 3258093488453501312L;

	private final RandomVariableFactory	abstractRandomVariableFactory;

	private final TimeDiscretization	simulationTimeDiscretization;
	private final TimeDiscretization	timeToMaturityDiscretization;

	private final Map<Integer, Map<Integer, Integer>> 	indexMap = new ConcurrentHashMap<>();
	private RandomVariable[] volatility;
	private final	boolean		isCalibrateable;

	public LIBORVolatilityModelPiecewiseConstantWithMercurioModification(final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization simulationTimeDiscretization, final TimeDiscretization timeToMaturityDiscretization, final RandomVariable[] volatility, final boolean isCalibrateable) {
		super(timeDiscretization, liborPeriodDiscretization);

		abstractRandomVariableFactory = new RandomVariableFromArrayFactory();

		/*
		 * Build index map
		 */
		

		final double maxMaturity = liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1);
		int volatilityIndex = 0;
		for(int simulationTime=0; simulationTime<simulationTimeDiscretization.getNumberOfTimes(); simulationTime++) {
			final HashMap<Integer, Integer> timeToMaturityIndexing = new HashMap<>();
			for(int timeToMaturity=0; timeToMaturity<timeToMaturityDiscretization.getNumberOfTimes(); timeToMaturity++) {
				if(simulationTimeDiscretization.getTime(simulationTime)+timeToMaturityDiscretization.getTime(timeToMaturity) > maxMaturity) {
					continue;
				}

				timeToMaturityIndexing.put(timeToMaturity,volatilityIndex++);
			}
			indexMap.put(simulationTime, timeToMaturityIndexing);
		}

		if(volatility.length == 1) {
			this.volatility = new RandomVariable[volatilityIndex];
			Arrays.fill(this.volatility, volatility[0]);
		}
		else if(volatility.length == volatilityIndex) {
			this.volatility = volatility.clone();
		}
		else {
			throw new IllegalArgumentException("Volatility length does not match number of free parameters.");
		}

		if(volatilityIndex != this.volatility.length) {
			throw new IllegalArgumentException("volatility.length should equal simulationTimeDiscretization.getNumberOfTimes()*timeToMaturityDiscretization.getNumberOfTimes().");
		}
		this.simulationTimeDiscretization = simulationTimeDiscretization;
		this.timeToMaturityDiscretization = timeToMaturityDiscretization;
		this.isCalibrateable = isCalibrateable;
	}

	public LIBORVolatilityModelPiecewiseConstantWithMercurioModification(final RandomVariableFactory abstractRandomVariableFactory, final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization simulationTimeDiscretization, final TimeDiscretization timeToMaturityDiscretization, final double[][] volatility, final boolean isCalibrateable) {
		super(timeDiscretization, liborPeriodDiscretization);

		this.abstractRandomVariableFactory = abstractRandomVariableFactory;

		/*
		 * Build index map
		 */
		final double maxMaturity = liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1);	
		int volatilityIndex = 0;
		for(int simulationTime=0; simulationTime<simulationTimeDiscretization.getNumberOfTimes(); simulationTime++) {
			final Map<Integer, Integer> timeToMaturityIndexing = new ConcurrentHashMap<>();
			for(int timeToMaturity=0; timeToMaturity<timeToMaturityDiscretization.getNumberOfTimes(); timeToMaturity++) {
				if(simulationTimeDiscretization.getTime(simulationTime)+timeToMaturityDiscretization.getTime(timeToMaturity) > maxMaturity) {
					continue;
				}

				timeToMaturityIndexing.put(timeToMaturity,volatilityIndex++);
			}
			indexMap.put(simulationTime, timeToMaturityIndexing);
		}

		// Flatten parameter matrix
		this.volatility = new RandomVariable[volatilityIndex];
		for(final Integer simulationTime : indexMap.keySet()) {
			for(final Integer timeToMaturity : indexMap.get(simulationTime).keySet()) {
				this.volatility[indexMap.get(simulationTime).get(timeToMaturity)] = abstractRandomVariableFactory.createRandomVariable(volatility[simulationTime][timeToMaturity]);
			}
		}

		this.simulationTimeDiscretization = simulationTimeDiscretization;
		this.timeToMaturityDiscretization = timeToMaturityDiscretization;
		this.isCalibrateable = isCalibrateable;
	}
	
	// ----------> IMPORTANTE
		// praticamente timeToMaturityDiscretization ti dice a quale LIBOR ci stiamo riferendeo quindi se prendi una discretizzazione annuale significa che, se i libor sono trimestrali, partendo dAL TEMPO 0, i primi 4 libor avranno lo stesso valore della volatiltà poi i secondi 4, avranno un'altra volatilità è cosi via.. rappresenta l'asse vertical della tua idea delle volatilità del LIBOR, mentre la simulationDiscretization è l'asse orizzontale, cioè ti dice ogni quanto i valori della volatilità cambiano, se ad esempio trimestrale, significa che ogni volatilità resta la stesa per un trimestre di simulazione.
		// ma per avere un piecewise constant trimestrale non basterebbe specificare quindi un timeToMaturityDiscretization? e optionMaturityDiscretization dovrebbe essere tipo (0,40).. NO! questo praticamente ti da una time-homogeneous piecewise constatn!!
		// NB:timeToMaturityDiscretization è proprio time to maturity quindi assumendo sia trimestrale un L(0,1;0) prenderà il valore della volatilità con Index=4.
		//  la reale piecewise constant volatility è data dalla matrice volatility[simulationTime][timeToMaturity])
	public LIBORVolatilityModelPiecewiseConstantWithMercurioModification(final RandomVariableFactory abstractRandomVariableFactory, final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization simulationTimeDiscretization, final TimeDiscretization timeToMaturityDiscretization, final double[] volatility, final boolean isCalibrateable) {
		super(timeDiscretization, liborPeriodDiscretization);

		this.abstractRandomVariableFactory = abstractRandomVariableFactory;

		/*
		 * Build index map 
		 */
		
		final double maxMaturity = liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1);
		int volatilityIndex = 0;
		for(int simulationTime=0; simulationTime<simulationTimeDiscretization.getNumberOfTimes(); simulationTime++) {
			final HashMap<Integer, Integer> timeToMaturityIndexing = new HashMap<>();
			for(int timeToMaturity=0; timeToMaturity<timeToMaturityDiscretization.getNumberOfTimes(); timeToMaturity++) {
				if(simulationTimeDiscretization.getTime(simulationTime)+timeToMaturityDiscretization.getTime(timeToMaturity) > maxMaturity) {
					continue;	
				}
				timeToMaturityIndexing.put(timeToMaturity,volatilityIndex++);
			}
			//quindi al simulation time 0 avrai un timeToMaturityIndexing che mi collega timeToMaturity 1..n a volatilityIndex 1..n, poi al simulation time 1 timeToMaturityIndexing  mi collega timeToMaturity 1..n a volatilityIndex n+1..2n
			indexMap.put(simulationTime, timeToMaturityIndexing);
			//pensa bene a come deve essere una piecewise volatility, è chiaro che ad ogni simulation time, il numero delle volatilità che devono essere salvate cambia. Ogni volta che sarà chiamato il getVolatility, ti ritornerà un vettore con i valori della volatilità, la cui dimensione è quella esatta, cioe solo per i LIBOR che non sono ancora fissati 
		}

		if(volatility.length == 1) {
			this.volatility = new RandomVariable[volatilityIndex];
			//assegna il valore di volatility[0] a tutto l'array
			Arrays.fill(this.volatility, abstractRandomVariableFactory.createRandomVariable(volatility[0]));
		}
		else if(volatility.length == volatilityIndex) {
			this.volatility = new RandomVariable[volatilityIndex];
			for(int i=0; i<volatility.length; i++) {
				this.volatility[i] = abstractRandomVariableFactory.createRandomVariable(volatility[i]);
			}
		}
		else {
			throw new IllegalArgumentException("Volatility length does not match number of free parameters.");
		}

		if(volatilityIndex != this.volatility.length) {
			throw new IllegalArgumentException("volatility.length should equal simulationTimeDiscretization.getNumberOfTimes()*timeToMaturityDiscretization.getNumberOfTimes().");
		}
		this.simulationTimeDiscretization = simulationTimeDiscretization;
		this.timeToMaturityDiscretization = timeToMaturityDiscretization;
		this.isCalibrateable = isCalibrateable;
	}

	public LIBORVolatilityModelPiecewiseConstantWithMercurioModification(final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization simulationTimeDiscretization, final TimeDiscretization timeToMaturityDiscretization, final double[] volatility, final boolean isCalibrateable) {
		this(new RandomVariableFromArrayFactory(), timeDiscretization, liborPeriodDiscretization, simulationTimeDiscretization, timeToMaturityDiscretization, volatility, isCalibrateable);
	}

	public LIBORVolatilityModelPiecewiseConstantWithMercurioModification(final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization simulationTimeDiscretization, final TimeDiscretization timeToMaturityDiscretization, final double volatility, final boolean isCalibrateable) {
		this(timeDiscretization, liborPeriodDiscretization, simulationTimeDiscretization, timeToMaturityDiscretization, new double[] { volatility }, isCalibrateable);
	}

	public LIBORVolatilityModelPiecewiseConstantWithMercurioModification(final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization simulationTimeDiscretization, final TimeDiscretization timeToMaturityDiscretization, final double[] volatility) {
		this(timeDiscretization, liborPeriodDiscretization, simulationTimeDiscretization, timeToMaturityDiscretization, volatility, true);
	}

	public LIBORVolatilityModelPiecewiseConstantWithMercurioModification(final TimeDiscretization timeDiscretization, final TimeDiscretization liborPeriodDiscretization, final TimeDiscretization simulationTimeDiscretization, final TimeDiscretization timeToMaturityDiscretization, final double volatility) {
		this(timeDiscretization, liborPeriodDiscretization, simulationTimeDiscretization, timeToMaturityDiscretization, new double[] { volatility });
	}

	@Override
	public RandomVariable[] getParameter() {
		if(isCalibrateable) {
			return volatility;
		} else {
			return null;
		}
	}

	@Override
	public LIBORVolatilityModel getCloneWithModifiedParameter(final RandomVariable[] parameter) {
		return new LIBORVolatilityModelPiecewiseConstantWithMercurioModification(
				super.getTimeDiscretization(),
				super.getLiborPeriodDiscretization(),
				simulationTimeDiscretization,
				timeToMaturityDiscretization,
				parameter,
				isCalibrateable
				);
	}

	@Override
	public RandomVariable getVolatility(final int timeIndex, final int liborIndex) {
		final double time             = getTimeDiscretization().getTime(timeIndex);	
//----->	
//MERCURIO:
//		final double maturity         = getLiborPeriodDiscretization().getTime(liborIndex);
		final double maturity         = getLiborPeriodDiscretization().getTime(liborIndex+1);
		final double timeToMaturity   = maturity-time;

		double volatilityInstanteaneous;
		if(timeToMaturity <= 0)
		{
			volatilityInstanteaneous = 0.0;   // This forward rate is already fixed, no volatility

			return abstractRandomVariableFactory == null ? new Scalar(0.0) : abstractRandomVariableFactory.createRandomVariable(time, volatilityInstanteaneous);
		}
		else
		{
			int timeIndexSimulationTime = simulationTimeDiscretization.getTimeIndexNearestLessOrEqual(time);
			if(timeIndexSimulationTime < 0) {
				timeIndexSimulationTime = 0;
			}
			if(timeIndexSimulationTime >= simulationTimeDiscretization.getNumberOfTimes()) {
				timeIndexSimulationTime--;
			}

			int timeIndexTimeToMaturity = timeToMaturityDiscretization.getTimeIndexNearestLessOrEqual(timeToMaturity);
			if(timeIndexTimeToMaturity < 0) {
				timeIndexTimeToMaturity = 0;
			}
			if(timeIndexTimeToMaturity >= timeToMaturityDiscretization.getNumberOfTimes()) {
				timeIndexTimeToMaturity--;
			}

			final int parameterIndex = indexMap.get(timeIndexSimulationTime).get(timeIndexTimeToMaturity);
			return volatility[parameterIndex];
		}
	}

	@Override
	public Object clone() {
		return new LIBORVolatilityModelPiecewiseConstantWithMercurioModification(
				super.getTimeDiscretization(),
				super.getLiborPeriodDiscretization(),
				simulationTimeDiscretization,
				timeToMaturityDiscretization,
				volatility.clone(),
				isCalibrateable
				);
	}


	/**
	 * @return the simulationTimeDiscretization
	 */
	public TimeDiscretization getSimulationTimeDiscretization() {
		return simulationTimeDiscretization;
	}

	/**
	 * @return the timeToMaturityDiscretization
	 */
	public TimeDiscretization getTimeToMaturityDiscretization() {
		return timeToMaturityDiscretization;
	}

	@Override
	public LIBORVolatilityModel getCloneWithModifiedData(final Map<String, Object> dataModified) {
		RandomVariableFactory abstractRandomVariableFactory = this.abstractRandomVariableFactory;
		TimeDiscretization timeDiscretization = this.getTimeDiscretization();
		TimeDiscretization liborPeriodDiscretization = this.getLiborPeriodDiscretization();
		TimeDiscretization simulationTimeDiscretization = this.getSimulationTimeDiscretization();
		TimeDiscretization timeToMaturityDiscretization = this.getTimeToMaturityDiscretization();
		double[][] volatility = new double[simulationTimeDiscretization.getNumberOfTimes()][timeToMaturityDiscretization.getNumberOfTimes()];
		for(int i = 0;i<volatility.length;i++) {
			for(int j = 0;j<volatility[i].length;j++) {
				volatility[i][j] = this.volatility[indexMap.get(i).get(j)].doubleValue();
			}
		}

		if(dataModified != null) {
			// Explicitly passed covarianceModel has priority
			abstractRandomVariableFactory = (RandomVariableFactory)dataModified.getOrDefault("randomVariableFactory", abstractRandomVariableFactory);
			timeDiscretization = (TimeDiscretization)dataModified.getOrDefault("timeDiscretization", timeDiscretization);
			liborPeriodDiscretization = (TimeDiscretization)dataModified.getOrDefault("liborPeriodDiscretization", liborPeriodDiscretization);
			simulationTimeDiscretization = (TimeDiscretization)dataModified.getOrDefault("simulationTimeDiscretization", simulationTimeDiscretization);
			timeToMaturityDiscretization = (TimeDiscretization)dataModified.getOrDefault("timeToMaturityDiscretization", timeToMaturityDiscretization);


			if(dataModified.getOrDefault("volatility", volatility) instanceof double[][]) {
				volatility = (double[][])dataModified.getOrDefault("volatility", volatility);
			}
			else {
				// TODO Implement handling for double[], double, RV[], RV
				throw new UnsupportedOperationException("volatility parameter type not supported.");
			}
		}

		final LIBORVolatilityModel newModel = new LIBORVolatilityModelPiecewiseConstantWithMercurioModification(abstractRandomVariableFactory, timeDiscretization, liborPeriodDiscretization, simulationTimeDiscretization, timeToMaturityDiscretization, volatility, isCalibrateable);
		return newModel;
	}

}
