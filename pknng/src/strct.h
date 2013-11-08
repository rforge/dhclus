typedef struct c_node{
	int id;			//ID number of vertex
	void *strt;		// here you put what ever you want....
	struct c_node *alpha;	//alpha and beta is for up and down or left ans rigth 
	struct c_node *beta;	//
} c_node;
