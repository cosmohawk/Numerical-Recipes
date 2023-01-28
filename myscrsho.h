template<class T>
void scrsho(T &fx) {
  /* Graph the function or functor fx over the prompted-for interval x1,x2.
     Query for another plot until the user signals satisfaction. */

  // Number of function evaluations for each plot.
  const Int RES=500;
  
  // Corners of plots in points.
  const Doub XLL=75., XUR=525., YLL=250., YUR=700.;

  char *plotfilename = tmpnam(NULL);
  VecDoub xx(RES), yy(RES);
  Doub x1, x2;
  Int i;

  for (;;) {
    Doub ymax = -9.99e99, ymin = 9.99e99, del;

    // Query for another plot, quit if x1=x2.
    cout << endl << "Enter x1 x2 (x1=x2 to stop):" << endl;
    cin >> x1 >> x2;
    if (x1==x2) break;

    // Evaluate the function at equal intervals.
    // Find the largest and smallest values.
    for (i=0; i<RES; i++) {
      xx[i] = x1 + i*(x2-x1)/(RES-1.);
      yy[i] = fx(xx[i]);
      if (yy[i] > ymax) ymax=yy[i];
      if (yy[i] < ymin) ymin=yy[i];
    }

    del = 0.05*((ymax-ymin)+(ymax==ymin ? abs(ymax) : 0.));

    /* Plot commands, following, are in PSplot syntax. You can substitute 
       commands for your favourite plotting package. */
    PSpage pg(plotfilename);
    PSplot plot(pg,XLL,XUR,YLL,YUR);
    plot.setlimits(x1,x2,ymin-del,ymax+del);
    plot.frame();
    plot.autoscales();
    plot.lineplot(xx,yy);
    if (ymax*ymin < 0.) plot.lineseg(x1,0.,x2,0.);
    plot.display();
  }
  remove(plotfilename);
}
