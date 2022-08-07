package edge;

import node.Node;

public class Edge {
	Node parent;
	Node child;
	int numberOfLoss;

	public Edge(Node parent, Node child) {

		this.parent = parent;
		this.child = child;
	}

	public Node getParent() {
		return parent;
	}

	public void setParent(Node parent) {
		this.parent = parent;
	}

	public Node getChild() {
		return child;
	}

	public void setChild(Node child) {
		this.child = child;
	}

	public int getNumberOfLoss() {
		return numberOfLoss;
	}

	public void setNumberOfLoss(int numberOfLoss) {
		this.numberOfLoss = numberOfLoss;
	}

	@Override
	public String toString() {
		return "(" + parent + "->" + child + ") numberOfLoss=" + numberOfLoss;
	}

	public void incrementLossCount() {
		numberOfLoss++;
	}

	public void decrementLossCount() {
		numberOfLoss--;

	}

}
