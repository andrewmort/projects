double delta_boundary(unsigned idx, int dimen, double dist) {
	double cost = 0.0;
	for(int i = 0; i < locations->size(); i++) {
		double xpos = locations->at(i).x;
		double ypos = locations->at(i).y;

		if(xpos < 0) cost += pow(xpos / alpha, 2);
		if(ypos < 0) cost += pow(ypos / alpha, 2);
		if(xpos > chipx) cost += pow((xpos - chipx) / alpha, 2);
		if(ypos > chipy) cost += pow((ypos - chipy) / alpha, 2);
	}
    	return cost;
}
