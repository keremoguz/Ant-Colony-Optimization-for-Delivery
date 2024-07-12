import java.util.ArrayList;
/**
 * The House class represents a house with its coordinates and properties.
 */
public class House {
    int id;
    double xCord;
    double yCord;
    ArrayList<House> neighbours = new ArrayList<House>();
    ArrayList<Double> distances = new ArrayList<Double>();
    ArrayList<Double> pheromoneLevels = new ArrayList<Double>();
    ArrayList<Double> edgeValues = new ArrayList<Double>();

    /**
     * Constructs a new house with default coordinates and ID.
     */
    House(){
        xCord = 0;
        yCord = 0;
        id = 1;
    }
    /**
     * Constructs a new house with specified coordinates and ID.
     *
     * @param xCord The x-coordinate of the house.
     * @param yCord The y-coordinate of the house.
     * @param id The unique identifier of the house.
     */
    House(double xCord, double yCord,int id){
        this.xCord = xCord;
        this.yCord = yCord;
        this.id = id;
    }
    public String toString() {
        return ""+id;
    }

}
