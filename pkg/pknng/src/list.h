struct lnode;

typedef struct lnode{
	unsigned long int itag;
	unsigned long int jtag;
	unsigned long int cnt;
	unsigned long int cnt_root;
	double dist;
	char mark;
	struct lnode* next;
	struct lnode* down;
	struct lnode* up;
} lnode;

// from extern inline void shuffle_index(int ,int *,sample *);
typedef struct{
	int idx;
	double val;
} sample;

extern inline int add_node(lnode*,unsigned long int,double);
extern inline void init_root_node(lnode*,unsigned long int);
extern inline lnode * add_root_node(lnode *,unsigned long int);
extern inline void add_ij_node(lnode*,unsigned long int,unsigned long int,double);
extern inline void merge(lnode *,lnode *);
