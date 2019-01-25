import java.util.Arrays;
import java.util.Random;

public class BioSystem {

    private int L, K, s, s_max;
    private double c, alpha, timeElapsed, dt, dx;
    private double D = 0.1; //diffusion constant

    private Microhabitat[] microhabitats;
    Random rand = new Random();
    private int initialPop = 100;

    public BioSystem(int L, int S, double alpha, double dt, double dx){

        this.L = L;
        this.s = S;
        this.s_max = S;
        this.alpha = alpha;
        this.dt = dt;
        this.dx = dx;
        this.microhabitats = new Microhabitat[L];
        this.timeElapsed = 0.;

        for(int i = 0; i < L; i++){
            double c_i = Math.exp(alpha*(double)i) - 1.;
            microhabitats[i] = new Microhabitat(S, c_i);
        }
        microhabitats[0].fillWithWildType(initialPop);
    }

    public int getL(){
        return L;
    }
    public int getS(){return s;}
    public double getTimeElapsed(){
        return timeElapsed;
    }

    public void updateMicrohabitats(Microhabitat[] updatedSystem){
        //this updates the system with the new values for S and N
        this.microhabitats = updatedSystem;
    }

    public Microhabitat getMicrohabitats(int i){
        return microhabitats[i];
    }

    public int getTotalN(){
        int runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN_alive()+m.getN_dead();
        }
        return runningTotal;
    }

    public int getCurrentLivePopulation(){
        int runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN_alive();
        }
        return runningTotal;
    }

    public int[] getLiveSpatialDistributionArray(){
        int[] mh_pops = new int[L];
        for(int i = 0; i < L; i++){
            mh_pops[i] = microhabitats[i].getN_alive();
        }
        return mh_pops;
    }

    public int[] getDeadSpatialDistributionArray(){
        int[] mh_pops = new int[L];
        for(int i = 0; i < L; i++){
            mh_pops[i] = microhabitats[i].getN_dead();
        }
        return mh_pops;
    }


    public double[] getGrowthRatesArray(){
        double[] mh_gRates = new double[L];
        for(int i = 0; i < L; i++){
            mh_gRates[i] = microhabitats[i].replication_or_death_rate();
        }
        return mh_gRates;
    }

    public int[] getNutrientsArray(){
        int[] mh_S = new int[L];
        for(int i = 0; i < L; i++){
            mh_S[i] = microhabitats[i].getS();
        }
        return mh_S;
    }


    public void performAction(){

        Microhabitat[] newSystem = microhabitats.clone();

        for(int i = 0; i < L; i++){

            int s_i = microhabitats[i].getS();
            int nalive_i = microhabitats[i].getN_alive();
            int ndead_i = microhabitats[i].getN_dead();

            if(microhabitats[i].itIsGrowing()){

                double growthRate_i = microhabitats[i].replication_rate();
                s_i += -dt*nalive_i*growthRate_i;


                if(i==0){
                    nalive_i += dt*(D*(microhabitats[i+1].getN_alive() + 0 - 2*nalive_i)/(dx*dx) + nalive_i*growthRate_i);
                }else if(i==L-1){
                    nalive_i += dt*(D*(0 + microhabitats[i-1].getN_alive() - 2*nalive_i)/(dx*dx)+ nalive_i*growthRate_i);
                }else{
                    nalive_i += dt*(D*(microhabitats[i+1].getN_alive() + microhabitats[i-1].getN_alive() - 2*nalive_i)/(dx*dx)
                            + nalive_i*growthRate_i);
                }

            }else{
                double deathRate_i = microhabitats[i].death_rate();
                //growth rate is -ve so no change in nutrients
                if(i==0){
                    double deltaN = dt*(D*(microhabitats[i+1].getN_alive() + 0 - 2*nalive_i)/(dx*dx) - nalive_i*deathRate_i);

                    nalive_i += deltaN;
                    ndead_i += Math.abs(deltaN);
                }else if(i==L-1){
                    double deltaN = dt*(D*(0 + microhabitats[i-1].getN_alive() - 2*nalive_i)/(dx*dx)- nalive_i*deathRate_i);

                    nalive_i += deltaN;
                    ndead_i += Math.abs(deltaN);
                }else{
                    double deltaN = dt*(D*(microhabitats[i+1].getN_alive() + microhabitats[i-1].getN_alive() - 2*nalive_i)/(dx*dx)
                            - nalive_i*deathRate_i);

                    nalive_i += deltaN;
                    ndead_i += Math.abs(deltaN);
                }
            }
            newSystem[i].setS(s_i);
            newSystem[i].setN_alive(nalive_i);
            newSystem[i].setN_dead(ndead_i);
        }


        updateMicrohabitats(newSystem);

        timeElapsed += dt;
    }


    public static void expGrad_popAndgRateDistbs(double biggestC, double dt, double dx){
        int L = 10, nReps = 2;
        int nMeasurements = 16;

        double duration = 10.;
        double interval = duration/nMeasurements;

        double alpha = BioSystem.calculateAlpha(L, biggestC);
        int S = 500;

        String filename_alive = "FGTA_death-alpha="+String.valueOf(alpha)+"-aliveDistribution-continuum";
        String filename_dead = "FGTA_death-alpha="+String.valueOf(alpha)+"-deadDistribution-continuum";
        String filename_gRate = "FGTA_death-alpha="+String.valueOf(alpha)+"-gRateDistribution-continuum";
    }


    public static void exponentialGradient_spatialAndGRateDistributions(double biggestC, double dt, double dx){

        int L = 10, nReps = 2;
        int nTimeMeasurements = 20;

        double duration = 2.;
        double interval = duration/(double)nTimeMeasurements;

        double alpha = BioSystem.calculateAlpha(L, biggestC);
        int S = 500;

        String filename_alive = "FGTA_death-alpha="+String.valueOf(alpha)+"-aliveDistribution-continuum";
        String filename_dead = "FGTA_death-alpha="+String.valueOf(alpha)+"-deadDistribution-continuum";
        String filename_gRate = "FGTA_death-alpha="+String.valueOf(alpha)+"-gRateDistribution-continuum";


        int[][][] allN_alive = new int[nReps][][];
        int[][][] allN_dead = new int[nReps][][];
        double[][][] allGRates = new double[nReps][][];

        for(int r = 0; r < nReps; r++){

            boolean alreadyRecorded = false;

            int[][] alivePopsOverTime = new int[nTimeMeasurements+1][];
            int[][] deadPopsOverTime = new int[nTimeMeasurements+1][];
            double[][] gRatesOverTime = new double[nTimeMeasurements+1][];
            int timerCounter = 0;


            BioSystem bs = new BioSystem(L, S, alpha, dt, dx);

            while(bs.timeElapsed <= duration+0.2*interval){

                if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.1*interval) && !alreadyRecorded){

                    System.out.println("rep: "+r+"\ttime elapsed: "+String.valueOf(bs.getTimeElapsed()));

                    //System.out.println(Arrays.toString(bs.getNutrientsArray()));

                    alivePopsOverTime[timerCounter] = bs.getLiveSpatialDistributionArray();
                    deadPopsOverTime[timerCounter] = bs.getDeadSpatialDistributionArray();
                    gRatesOverTime[timerCounter] = bs.getGrowthRatesArray();

                    alreadyRecorded = true;
                    timerCounter++;
                    System.out.println(timerCounter);
                }
                if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;


                bs.performAction();
            }

            //System.out.println(Arrays.deepToString(alivePopsOverTime));
            allN_alive[r] = alivePopsOverTime;
            allN_dead[r] = deadPopsOverTime;
            allGRates[r] = gRatesOverTime;



        }
        //System.out.println(Arrays.toString(allN_alive));
        double[][] averagedAlivePopDistributions = Toolbox.averagedResults(allN_alive);
        double[][] averagedDeadPopDistributions = Toolbox.averagedResults(allN_dead);
        double[][] averagedGRateDistributions = Toolbox.averagedResults(allGRates);

        //print the live and dead results to two separate files, then just join them together in
        //gnuplot or something
        Toolbox.printAveragedResultsToFile(filename_alive, averagedAlivePopDistributions);
        Toolbox.printAveragedResultsToFile(filename_dead, averagedDeadPopDistributions);
        Toolbox.printAveragedResultsToFile(filename_gRate, averagedGRateDistributions);

    }


    public static double calculateAlpha(int L, double c_max){

        return Math.log(c_max+1)/L;
    }

}
