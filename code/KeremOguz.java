/**
 * This class represents a solution to the Traveling Salesman Problem (TSP) using either a brute-force method or Ant Colony Optimization (ACO).
 * It aims to find the shortest path that visits a given set of houses exactly once and returns to the starting point.
 * The class includes methods for both brute-force and ACO algorithms, as well as utility methods for path calculation,
 * pheromone update, edge value calculation, and permutation generation.
 *
 * @author Kerem Oguz ID : 2022400270
 * @since 11.05.2024
 */


import javax.print.attribute.HashPrintJobAttributeSet;
import java.awt.*;
import java.io.DataInput;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;

public class KeremOguz {
    public static double minimumDistance = Double.MAX_VALUE; // minimum distance of the path is kept
    public static ArrayList<House> optimalPath = new ArrayList<>(); // path with minimum distant is kept

    public static void main(String[] args) throws FileNotFoundException {
        long startTime = System.currentTimeMillis();
        int chosenMethod = 2;// 1 for brute-force, 2 for Ant Colony Optimization
        int whichToDraw = 2; // 1 for Shortest Path, 2 for Pheromone Intensity Map
        int iterationCount = 100;
        int antCountPerIteration = 50;
        double degradationFactor = 0.9;
        double alpha = 1.8; // alpha and beta values are chosen accordingly for an accurate solution
        double beta = 1.5;
        double pheromoneIntensity = 0.1;
        double qValue = 0.0001;

        ArrayList<House> houses = new ArrayList<>(); // houses including Migros is kept
        String houseCoords = "input05.txt"; // the name of the input file with house coordinates
        File houseCoordsFile = new File(houseCoords);
        if (!houseCoordsFile.exists()) {
            System.out.printf("%s cannot be found.", houseCoordsFile);
            System.exit(1);
        }

        Scanner houseCoordsInputFile = new Scanner(houseCoordsFile);

        int id = 0; // id denotes the order of the house (for Migros it is 1)
        while (houseCoordsInputFile.hasNextLine()) {
            id++;// for every new input line it gets incremented
            String line = houseCoordsInputFile.nextLine(); // current line is kept
            String[] info = line.split(","); // since line is separated with commas, we parse it
            double xCord = Double.parseDouble(info[0].trim()); // the first element is the x coordinate of the house
            double yCord = Double.parseDouble(info[1].trim()); // the second element is the y coordinate of the house
            houses.add(new House(xCord, yCord,id)); // new house object is constructed with the x,y coordinates and id data fields
        }
        houseCoordsInputFile.close(); // for safety the file is closed

        for (House houseI : houses) {
            houseI.neighbours = new ArrayList<>(houses); // ith house's neighbours are created
            // (house itself is also considered as a neighbour with distance 0 to itself)
            for (House houseJ : houses) {
                double distance = distance(houseI, houseJ); // distance between the house and the neighbour is calculated
                houseI.distances.add(distance);
                houseI.pheromoneLevels.add(pheromoneIntensity); // added initial pheromone intensity
                houseI.edgeValues.add(edgeValCalculator(distance, pheromoneIntensity, alpha, beta));// calculated edge value is added to edge values of the given house
            }
        }


        if(chosenMethod==1) {
            // if the method chosen is 1 brute-force method will be applied
            int[] arr = new int[houses.size()-1]; // array that is going to be permuted is initialized with size-1
            // first house(which is actually Migros) is not included in the permutation array to reduce complexity from n! to (n-1)!
            // since they constitute a circle not including the source does not change the result
            for (int i = 0; i<arr.length;i++){
                arr[i] = i+1; // excluding 0th element(Migros)
            }
            bruteForce(houses,arr,0);
            optimalPath.add(houses.getFirst()); // adding the Migros as the first and last element to close the circle
            optimalPath.add(0,houses.getFirst());

        } else if (chosenMethod == 2) {
            // if the method chosen is 2 Ant Colony Optimization will be applied
            // the method below returns the optimized path
            optimalPath = antColonyOptimization(houses, iterationCount, antCountPerIteration, qValue, degradationFactor, alpha, beta);
            // since the source is arbitrary chosen the resulting array will start with a random house,
            // to reconstruct the same path but starting with Migros, before and after Migros array lists are created
            ArrayList<House>beforeMigros = new ArrayList<>();
            ArrayList<House>afterMigros = new ArrayList<>();
            boolean isAfterMigros = false;
            for(House house: optimalPath){
                if(house.toString().equals("1")) isAfterMigros = true; // meaning we are at Migros so,from now add house to afterMigros
                if(isAfterMigros) afterMigros.add(house);
                else beforeMigros.add(house);
            }
            for(int i = 1; i< beforeMigros.size();i++){
                //adds the elements before the Migros to after the migros which reconstructs the path starting with Migros
                afterMigros.add(beforeMigros.get(i));
            }
            afterMigros.add(afterMigros.getFirst()); // closing the circle by adding Migros as last element
            optimalPath = afterMigros; // reassigning the reconstructed path to optimal path
        }
        else {
            System.out.println("No valid chosen method value."); // if neither 1 nor 2 is not chosen, informs the user then exist the programme
            System.exit(0);
        }


        // Visualization
        StdDraw.enableDoubleBuffering(); // it is used for smooth graphing
        StdDraw.setCanvasSize(700, 700);
        StdDraw.setXscale(0, 1); // x and y axis are scaled
        StdDraw.setYscale(0, 1);
        StdDraw.setPenColor(StdDraw.BLACK);
        if(whichToDraw == 2) {
            // draws pheromone intensity map
            for (int i = 1; i < houses.size(); i++) {
                House house = houses.get(i);
                for (int j = 0; j < house.neighbours.size(); j++) {
                    double pheromoneLevel = house.pheromoneLevels.get(j);
                    StdDraw.setPenRadius(pheromoneLevel*0.88);// pen thickness is adjusted with respect to pheromone level
                    StdDraw.line(house.xCord, house.yCord, house.neighbours.get(j).xCord, house.neighbours.get(j).yCord); // line is drawn between two houses
                }
            }
        }
        if (whichToDraw == 1) {
            StdDraw.setPenColor(StdDraw.BLACK);
            StdDraw.setPenRadius(0.005);
            // draws the optimal path only
            for (int i = 0; i < optimalPath.size() - 1; i++) {
                StdDraw.line(optimalPath.get(i).xCord, optimalPath.get(i).yCord, optimalPath.get(i + 1).xCord, optimalPath.get(i + 1).yCord);
            }
        }
        Font font = new Font("Arial", Font.PLAIN, 10); // font is set to fit in the circles
        StdDraw.setFont(font);
        for (House house : houses) {
            // if id = 1 (i.e Migros) denote it by drawing in orange, otherwise draw it to gray
            if(house.id == 1 && whichToDraw == 1) StdDraw.setPenColor(StdDraw.PRINCETON_ORANGE);
            else StdDraw.setPenColor(StdDraw.LIGHT_GRAY);

            StdDraw.filledCircle(house.xCord, house.yCord, 0.017);
            StdDraw.setPenColor(StdDraw.BLACK);
            StdDraw.text(house.xCord, house.yCord, house.toString());
        }
        StdDraw.show(); // show the graph on the canvas

        // informing the user with the method name, runtime, optimal distance and the path itself

        if (chosenMethod == 1) System.out.println("Method: Brute-Force Method");
        else System.out.println("Method : Ant Colony Optimization");
        System.out.printf("Shortest Distance: %.5f \n",pathDistance(optimalPath));
        System.out.println("Shortest Path: "+optimalPath);
        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        String secDuration =  duration/1000+ "."+ duration%1000;
        System.out.println("Time it takes to find the shortest path: " + secDuration + " seconds.");
    }

    /**
     Performs Ant Colony Optimization to find the shortest path among a given set of houses.

     @param houses The list of houses among which the shortest path is to be found.
     @param iterationCount The number of iterations to perform the optimization.
     @param antCountPerIteration The number of ants (traversals) to be simulated in each iteration.
     @param qValue Q-value representing the intensity of pheromone updates.
     @param degradationFactor Rate at which pheromones degrade after each iteration.
     @param alpha Parameter influencing the importance of pheromone levels in path selection.
     @param beta Parameter influencing the importance of heuristic information in path selection.
     @return The list of houses representing the shortest path found by the algorithm.
     */
    public static ArrayList<House> antColonyOptimization(ArrayList<House> houses, int iterationCount, int antCountPerIteration, double qValue, double degradationFactor, double alpha, double beta) {
        double minDist = Double.MAX_VALUE;
        ArrayList<House> minDistPath = new ArrayList<>();
        Random rand = new Random();
        for (int i = 0; i < iterationCount; i++) {
            for (int j = 0; j < antCountPerIteration; j++) {
                int randomIndex = rand.nextInt(houses.size());
                House source = houses.get(randomIndex);
                ArrayList<House> currentPath = antTraverse(source);
                double pathDist = pathDistance(currentPath);
                updatePheromone(currentPath, pathDist, qValue, alpha, beta);
                if (minDist > pathDist) {
                    minDist = pathDist;
                    minDistPath = new ArrayList<>(currentPath);
                }
            }
            degradePheromones(houses, degradationFactor);
        }
        return minDistPath;
    }

    /**
     * Conducts an ant traverse algorithm on the given source house.
     * The algorithm randomly selects adjacent houses to traverse based on edge values.
     *
     * @param source The starting house for the traversal.
     * @return An ArrayList containing the houses visited during the traversal,
     *         including the starting house at both the beginning and end.
     */
    public static ArrayList<House> antTraverse(House source) {
        Random rand = new Random();
        ArrayList<House> visited = new ArrayList<>(); // visited houses are kept
        while (true) {
            visited.add(source);
            double totalEdgeVal = calcTotalEdgeValue(source, visited); // recalculate the total edge value for the current house
            if (visited.size() == source.neighbours.size()) {
                // if all houses are visited, break
                break;
            }
            double[] probabilities = new double[source.neighbours.size()]; // probilities of the houses are kept
            double sum = 0;
            for (int i = 0; i < source.neighbours.size(); i++) {
                if (!visited.contains(source.neighbours.get(i))) {
                    probabilities[i] = source.edgeValues.get(i) / totalEdgeVal;
                    sum += probabilities[i]; // total probibility is accumulated
                }
            }
            double randNum = rand.nextDouble() * sum; // random number between 0 and total possibilty is created
            sum = 0;
            for (int i = 0; i < source.neighbours.size(); i++) {
                if (!visited.contains(source.neighbours.get(i))) {
                    sum += probabilities[i];
                    if (randNum <= sum) {
                        // if the random number lies in the interval visit there
                        source = source.neighbours.get(i);
                        break;
                    }
                }
            }
        }
        visited.add(visited.get(0)); // Return to the starting point
        return visited;
    }

    /**
     * Updates the total edge value while houses are being visited
     *
     * @param source The node being investigated
     * @param visited Visited houses so far
     * @return Updated total edge value for that node
     */
    public static double calcTotalEdgeValue(House source, ArrayList<House> visited) {
        double initialEdgeValue = sum(source.edgeValues);
        for (House house : visited) {
            initialEdgeValue -= house.edgeValues.get(house.neighbours.indexOf(source));
        }
        return initialEdgeValue;
    }

    /**
     * Calculates the distance of an unclosed path
     *
     * @param path The path consisting of houses
     * @return The total distance
     */
    public static double pathDistance(ArrayList<House> path) {
        double pathDist = 0;
        for (int i = 0; i < path.size() - 1; i++) {
            pathDist += distance(path.get(i), path.get(i + 1));
        }
        return pathDist;
    }

    /**
     * Updates the pheromone levels on the edges of the houses based on the current ant path.
     * Also updates the edge values using an edge value calculator function.
     *
     * @param currentPath An ArrayList containing the houses visited by the ant during its traversal.
     * @param pathDist The distance traveled by the ant along the path.
     * @param qValue The pheromone quantity constant.
     * @param alpha The alpha parameter used in the edge value calculation.
     * @param beta The beta parameter used in the edge value calculation.
     */
    public static void updatePheromone(ArrayList<House> currentPath, double pathDist, double qValue, double alpha, double beta) {
        double delta = qValue / pathDist;
        for (int i = 0; i < currentPath.size() - 1; i++) {
            House firstHouse = currentPath.get(i);
            House secondHouse = currentPath.get(i + 1);
            int index1 = firstHouse.neighbours.indexOf(secondHouse);
            int index2 = secondHouse.neighbours.indexOf(firstHouse);
            double prevVal1 = firstHouse.pheromoneLevels.get(index1);
            double prevVal2 = secondHouse.pheromoneLevels.get(index2);

            // updated the pheremone levels by adding delta
            firstHouse.pheromoneLevels.set(index1, prevVal1 + delta);
            secondHouse.pheromoneLevels.set(index2, prevVal2 + delta);

            // since pheromone levels are updated, edge values are updated accordingly
            double d = distance(firstHouse, secondHouse);
            double lastph = firstHouse.pheromoneLevels.get(index1);
            double edgeVal = edgeValCalculator(d, lastph, alpha, beta);
            firstHouse.edgeValues.set(index1, edgeVal);
            secondHouse.edgeValues.set(index2, edgeVal);
        }
    }
    /**
     * Degrades the pheromone levels on the edges of all houses by a given degradation factor.
     *
     * @param houses An ArrayList containing all the houses whose pheromone levels need to be degraded.
     * @param degradationFactor The factor by which the pheromone levels are degraded.
     */
    public static void degradePheromones(ArrayList<House> houses, double degradationFactor) {
        for (House house : houses) {
            for (int i = 0; i < house.pheromoneLevels.size(); i++) {
                double prevPheromone = house.pheromoneLevels.get(i);
                house.pheromoneLevels.set(i, prevPheromone * degradationFactor);
            }
        }
    }

    /**
     * The distance between two house objects is calculated
     *
     * @param house1 House 1 object
     * @param house2 House 2 object
     * @return the distance
     */
    public static double distance(House house1, House house2) {
        return Math.sqrt(Math.pow((house1.xCord - house2.xCord), 2) + Math.pow((house1.yCord - house2.yCord), 2));
    }
    /**
     * Calculates the edge value using the provided distance, pheromone intensity, alpha, and beta parameters.
     *
     * @param distance The distance between the houses connected by the edge.
     * @param pheromoneIntensity The intensity of pheromone on the edge.
     * @param alpha The alpha parameter used in the calculation.
     * @param beta The beta parameter used in the calculation.
     * @return The calculated edge value.
     */
    public static double edgeValCalculator(double distance, double pheromoneIntensity, double alpha, double beta) {
        if (distance == 0) return 0;
        return Math.pow(pheromoneIntensity, alpha) / Math.pow(distance, beta);
    }

    /**
     * It simply calculates the sum of Array list containing Doubles
     * @param arrayList ArrayList with Double values
     * @return sum of the ArrayList
     */
    public static double sum(ArrayList<Double> arrayList) {
        double total = 0;
        for (Double num : arrayList) {
            total += num;
        }
        return total;
    }

    /**
     * Generates permutations of house indices using a brute force approach to find the minimum path.
     *
     * @param houses An ArrayList containing all the houses to be visited.
     * @param arr An array representing the indices of houses.
     * @param k The current index in the permutation process.
     */
    private static void bruteForce(ArrayList<House> houses, int[] arr, int k) {
        if (k == arr.length) {
            checkIfMinPath(houses, arr);
        } else {
            for (int i = k; i < arr.length; i++) {
                int temp = arr[i];
                arr[i] = arr[k];
                arr[k] = temp;
                bruteForce(houses, arr, k + 1);
                temp = arr[k];
                arr[k] = arr[i];
                arr[i] = temp;
            }
        }
    }
    /**
     * Checks if the current path, represented by the provided array of house indices,
     * is shorter than the minimum recorded distance. If so, updates the minimum distance
     * and optimal path accordingly.
     *
     * @param houses An ArrayList containing all the houses to be visited.
     * @param arr An array representing the indices of houses in the current path.
     */
    public static void checkIfMinPath(ArrayList<House> houses, int[] arr) {
        double dist = 0;
        ArrayList<House> tempHouses = new ArrayList<>();
        for (int i = 0; i < arr.length - 1; i++) {
            House house1 = houses.get(arr[i]);
            House house2 = houses.get(arr[i + 1]);
            dist += distance(house1, house2);
        }
        dist+= distance(houses.getFirst(),houses.get(arr[0])); // distance to source (Migros) is also considered
        dist+= distance(houses.getFirst(),houses.get(arr[arr.length-1]));

        if (dist < minimumDistance) {
            // updated minimum distance and the path if the calculated distance is smaller then the previous min
            for(int i = 0 ; i<arr.length; i++){
                tempHouses.add(houses.get(arr[i]));
            }
            minimumDistance = dist;
            optimalPath = tempHouses;
        }
    }
}

