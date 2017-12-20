 John Hanson, 10/25/2017
Program "DepletionExtrapolation1"

This program is designed to estimate kcat/Km from a single depletion kinetic run.
In the simplest mode, if you choose 1 cycle, it will fit the data, assuming 
first order behavior and output a plot of the original data with the fit
superimposed, along with the parameters associated with the fit. The program
asks for the enzyme concentration, in molar, in order to calculate the kcat/Km
from the k(obs). If you don't have the enzyme concentration, you can just enter
1 and the kcat/Km will correspond to the k(obs).

But beware, the calculated kcat/Km is only accurate under conditions 
where [S] << Km (e.g. [S]/Km < 0.1). However, it is possible to determine the 
kcat/Km at various concentrations and then extrapolate the kcat/Km to [S] = 0.
See Crompton, I.E., and Waley, S.G. Biochem. J. "The determination of 
specificity contants in enzyme-catalyzed reactions" Biochem. J. 1986, 239,
221-224. However, it should also be possible to extrapolate to [S] = 0, by
using a single run -- determining kcat/Km at different starting points along
the reaction progress curve corresponds to depletion experiments at different [S].
Thus, for example, you should be able to determine kcat/Km for the full run,
then for the last 90% of the run, then the last 80% etc., and in this way
get a series of kcat/Km at different [S] that can then be used to extrapolate 
to the true kcat/Km.

The data should be in a tab delimited text file with the data in pairs of time 
(in seconds) followed by the absorbance (or any other variable that is proportional
to the amount of product. No other data (e.g., column headings) should be 
present. A convenient way to get the data in the appropriate form is to take
the data as two columns in an Excel file (most instrunments allow you to save
data in Excel format, and save as a tab delimited text file.

The program uses two different equations to fit the reaction progress curve:

     Abs = Abs(initial) + (1-exp(-k(obs)t))deltaAbs (This is actually Method 2)

Where deltaAbs is the Abs(final)-Abs(initial)
This to me is the conceptually simplest approach

However, another mathematically equivalent description of a first order reaction
progress curve is:

       Abs = Abs(0) + v(initial)*(1-exp(-k(obs)t))/k(obs) (This is Method 1)

 Where v(initial) = the initial rate

Both of these approaches should in theory give the same value for k(obs). And
in my experience they do. But I have left them both in for completeness, and
also because the second formulation is easily modified to fit to first order
processes that don't have a final rate of zero (e.g., burst kinetics).
INSERT CITATION and EQUATION

After k(obs) is determined, it is trivial to calculate an estimate of 
kcat/Km since k(obs)/[E] = kcat/Km. But just to reiterate, the reaction 
progress curve for an enzyme following simple Michaelis-Menten kinetics will
only be first order when [S] << Km. Thus if you do the experiment with [S]=Km,
you will see that the program does not provide a good fit. Nevertheless, as 
shown in the Crompton and Waley paper cited above, a plot of 1/k(obs) (or Km/kcat) 
versus [S] should give a reasonably linear relationship that can be used to 
extrapolate to [S] = 0.

This program does not ask for the [S] or Km, and in fact it is not necessary to
know these (nor the delta epsilon). The program plots kcat/Km versus relative [S], 
where the initial concentration of [S] is arbitrarily set at 1.
This expedient has the benefit of making the slope/intercept a reflection
of where [S] was relative to Km initially. 


* Gitter chatroom for ups_enzyme_kinetics development:

[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ups_enzyme_kinetics?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
https://gitter.im/ups_enzyme_kinetics/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link

[![Join the chat at https://gitter.im/phaustin/ups_enzyme_kinetics
](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/phaustin/?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
