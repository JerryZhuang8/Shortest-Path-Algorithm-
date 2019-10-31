import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

public class Graph {

  // Keep a fast index to nodes in the map
  private Map<Integer, Vertex> vertexNames;

  /**
   * Construct an empty Graph with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Graph() {
    vertexNames = new HashMap<>();
  }

  /**
   * Adds a vertex to the graph. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the graph
   */
  public void addVertex(Vertex v) {
    if (vertexNames.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertexNames.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the graph
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the graph
   */
  public Collection<Vertex> getVertices() {
    return vertexNames.values();
  }

  /**
   * Gets the vertex object with the given name
   * 
   * @param name
   *          (String) name of the vertex object requested
   * @return (Vertex) vertex object associated with the name
   */
  public Vertex getVertex(String name) {
    return vertexNames.get(Integer.parseInt(name));
  }

  /**
   * Adds a directed edge from vertex u to vertex v
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addEdge(int nameU, int nameV, Double cost) {
    if (!vertexNames.containsKey(nameU))
      throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
    if (!vertexNames.containsKey(nameV))
      throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
    Vertex sourceVertex = vertexNames.get(nameU);
    Vertex targetVertex = vertexNames.get(nameV);
    Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
    sourceVertex.addEdge(newEdge);
  }

  /**
   * Adds an undirected edge between vertex u and vertex v by adding a directed
   * edge from u to v, then a directed edge from v to u
   * 
   * @param name
   *          (String) name of vertex u
   * @param name2
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(int name, int name2, double cost) {
    addEdge(name, name2, cost);
    addEdge(name2, name, cost);
  }


  /**
   * Computes the euclidean distance between two points as described by their
   * coordinates
   * 
   * @param ux
   *          (double) x coordinate of point u
   * @param uy
   *          (double) y coordinate of point u
   * @param vx
   *          (double) x coordinate of point v
   * @param vy
   *          (double) y coordinate of point v
   * @return (double) distance between the two points
   */
  public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
    return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2));
  }

  /**
   * Computes euclidean distance between two vertices as described by their
   * coordinates
   * 
   * @param u
   *          (Vertex) vertex u
   * @param v
   *          (Vertex) vertex v
   * @return (double) distance between two vertices
   */
  public double computeEuclideanDistance(Vertex u, Vertex v) {
    return computeEuclideanDistance(u.x, u.y, v.x, v.y);
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanCost method.
   */
  public void computeAllEuclideanDistances() {
    for (Vertex u : getVertices())
      for (Edge uv : u.adjacentEdges) {
        Vertex v = uv.target;
        uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y);
      }
  }



  // STUDENT CODE STARTS HERE

  public void generateRandomVertices(int n) {
    vertexNames = new HashMap<>(); // reset the vertex hashmap
    Vertex vertexTemp = null; 
    for(int index = 0; index < n; index++) { 
         
         vertexTemp = new Vertex(index, (int)(101 * Math.random()), (int)(101 * Math.random())); 
         addVertex(vertexTemp); 
         
    
    }
    // Your code here...
    // now, we will add an edge so that you can go to and from any node to another node. 
    // 
    ArrayList<Vertex> modifier = new ArrayList<Vertex>(getVertices()); 
    for(int index = 0; index < n; index++) { 
        for(int innerIndex = index+1; innerIndex < n; innerIndex++) { //prevents repeats. kind of like a choose 2 algorithm. 
            addUndirectedEdge(index, innerIndex, 0.0);
        }
    } 
    computeAllEuclideanDistances(); // compute distances
  }

  public List<Edge> nearestNeighborTsp() {
     //randomize the beginning start point. 
     Vertex[] allVertices = getVertices().toArray(new Vertex[getVertices().size()]); 
     ArrayList<Vertex> notVisited = new ArrayList<Vertex>(); 
     for(int indexFiller = 0; indexFiller < allVertices.length; indexFiller++) { 
         notVisited.add(allVertices[indexFiller]); 
     }
     int scaler = allVertices.length; 
     Vertex cityOrigin = getVertex(Integer.toString((int)(Math.random()*scaler)));
     System.out.println("OKAY! LET US START AT VERTEX " + cityOrigin.name + ".");
     Vertex saver = cityOrigin; 
     notVisited.remove(saver); 
     Vertex tempClosest = null; 
     ArrayList<Edge> path = new ArrayList<Edge>();
     double minDistance = Double.MAX_VALUE; 
     int count = 0; 
     while(!notVisited.isEmpty()) { 
     for(int index = 0; index < cityOrigin.adjacentEdges.size(); index++) { 
         if(notVisited.contains(cityOrigin.adjacentEdges.get(index).target)) { 
            if(cityOrigin.adjacentEdges.get(index).distance < minDistance || count == 0) { 
               tempClosest = cityOrigin.adjacentEdges.get(index).target;
               minDistance = cityOrigin.adjacentEdges.get(index).distance; 
               count++; 
            }
         
         }
     
     
     }
     path.add(new Edge(cityOrigin, tempClosest, minDistance)); 
     notVisited.remove(tempClosest); 
     cityOrigin = tempClosest; 
     count = 0; 
     }
    path.add(new Edge(saver, cityOrigin, computeEuclideanDistance(saver, cityOrigin)));
    return path; // replace this line
  }

  public List<Edge> bruteForceTsp() {
    String permuteVertex = ""; 
    Vertex[] allVertices = getVertices().toArray(new Vertex[getVertices().size()]); 
    int graphSize = allVertices.length; 
    double tempLengthPath = 0; 
    ArrayList<Edge> idealPath = new ArrayList<Edge>(); 
    double minLengthPath = Double.MAX_VALUE; 
    for(int index = 0; index < allVertices.length; index++) { 
       permuteVertex = permuteVertex + " " + allVertices[index].name; 
    }
    String allCombos = combinations(permuteVertex);
   //returns all permutations. 
    int tempCount = 0; 
    ArrayList<Vertex> path = new ArrayList<Vertex>(); 
    for(int counter = 0; counter < allCombos.length(); counter++) { 
        if(allCombos.charAt(counter)==' ' || counter==allCombos.length()-1) { 
            tempCount++; 
            if(counter == allCombos.length()-1) { 
                path.add(getVertex(allCombos.substring(allCombos.length()-1, allCombos.length()))); 
            }
            if(tempCount > graphSize || counter == allCombos.length()-1) { 
            for(int n = 0; n < path.size(); n++)
            for(int checkTotal = 0; checkTotal < path.size()-1; checkTotal++) { 
                tempLengthPath += computeEuclideanDistance(path.get(checkTotal), path.get(checkTotal+1)); 
           }
            tempLengthPath += computeEuclideanDistance(path.get(0), path.get(path.size()-1));
            if(tempLengthPath < minLengthPath) { 
                minLengthPath = tempLengthPath; 
                idealPath.clear(); 
                for(int innerIndex = 0; innerIndex < path.size()-1; innerIndex++){
                idealPath.add(new Edge(path.get(innerIndex), path.get(innerIndex+1), computeEuclideanDistance(path.get(innerIndex),path.get(innerIndex+1)))); } 
                idealPath.add(new Edge(path.get(path.size()-1),path.get(0),computeEuclideanDistance(path.get(0),path.get(path.size()-1))));
                
            }
            path.clear();
            tempLengthPath = 0; 
            //check is done, we move on to the next permutation. 
            tempCount = 1; 
            }
        }
        else if(allCombos.charAt(counter)!=' ' && tempCount <= graphSize) { 
        for(int inner = counter; inner < allCombos.length(); inner++ ) { 
            if(allCombos.charAt(inner) == ' ' || inner == allCombos.length()-1) { 
                if(inner == allCombos.length()-1) { 
                path.add(getVertex((allCombos.substring(counter, inner+1)))); 
                }
                else { 
                path.add(getVertex(allCombos.substring(counter,inner)));
                }
                counter = inner-1; //so we start counting the next space again. 
                break; 
            }
        
        }      
            }
        
    }
     
    return idealPath; 
  
  }
  
  private String combinations(String combination) { 
      return permutation(combination, new String()); 
  
  
  }
  private String permutation(String permute, String finalString) { //string of all of the names. 
       //base case: determine if it is size one. the entered string will be of form "space1space2space3..." where 1,2,3 represent vertex names. 
       int numbVertices = 0; 
       for(int index = 0; index < permute.length(); index++) { 
           if(permute.charAt(index) != ' ') { 
              for(int counter = index; counter < permute.length(); counter++) { 
                   if(permute.charAt(counter) == ' ' || counter == permute.length()-1) { 
                       index = counter; 
                       numbVertices++; 
                       break; 
                   }
              
              }
           
           }
       
       }
       if(numbVertices == 1) { 
            return permute;  //BASE CASE. 
       }
       else { 
           String oneRemoved; 
           String tempString = ""; 
           String idealString = "";
           int spaceCounter = 0; 
           int count = 0; 
           String removedElement = ""; 
           for(int index = 0; index < permute.length(); index++) { 
           
                 if(permute.charAt(index)== ' ') { 
                 oneRemoved = new String(permute); //needed because we want to preserve the original String permute. 
                 for(int innerIndex = index+1; innerIndex < oneRemoved.length(); innerIndex++) { 
                      if(oneRemoved.charAt(innerIndex) == ' ' || innerIndex == oneRemoved.length()-1) { 
                         if(innerIndex == oneRemoved.length()-1) { 
                             removedElement = oneRemoved.substring(index+1, innerIndex+1); 
                             oneRemoved = oneRemoved.substring(0, index) + oneRemoved.substring(innerIndex+1, oneRemoved.length());
                         }
                         else if(oneRemoved.charAt(innerIndex)==' ') { 
                             removedElement = oneRemoved.substring(index+1, innerIndex); //exclude the space.
                             oneRemoved = oneRemoved.substring(0, index) + oneRemoved.substring(innerIndex, oneRemoved.length()); 
                             //once again, exclude the space here. 
                         }
                         break; 
                      } //here, end of the process of removing one element. 
                       
                 }
                 spaceCounter = -1; 
                 count = numbVertices-1; 
                  //Goal: Now, we create a string of all the permutations that begin with that removed element. 
                  //assign tempString to the recursive call that gives all the permutations on oneRemoved. 
                 tempString = permutation(oneRemoved, new String()); 
                 //now, it is necessary to add the removed element back to each possible permutation 
                 for(int insertIndex = 0; insertIndex < tempString.length(); insertIndex++) { 
                     if(tempString.charAt(insertIndex)== ' ') {  //space detected
                         spaceCounter++; 
                         if(spaceCounter % count == 0) { 
                              idealString = idealString + " " + removedElement + " ";                                                 
                         }
                         else { 
                             idealString = idealString + " ";
                         }                       
                     }
                     else { 
                        idealString = idealString + Character.toString(tempString.charAt(insertIndex)); 
                     
                     }
                               
                 }
                 finalString = finalString + idealString; 
                 idealString = ""; //necessary to reset for the next call. 
                 
                 
                 }
                
           }
           return finalString; 
          
       }
  
  
  
  
  
  }

  // STUDENT CODE ENDS HERE



  /**
   * Prints out the adjacency list of the graph for debugging
   */
  public void printAdjacencyList() {
    for (int u : vertexNames.keySet()) {
      StringBuilder sb = new StringBuilder();
      sb.append(u);
      sb.append(" -> [ ");
      for (Edge e : vertexNames.get(u).adjacentEdges) {
        sb.append(e.target.name);
        sb.append("(");
        sb.append(e.distance);
        sb.append(") ");
      }
      sb.append("]");
      System.out.println(sb.toString());
    }
  }
}
