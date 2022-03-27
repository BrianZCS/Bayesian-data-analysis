# Bayesian-data-analysis
Caleb Stevens, Rauf Salamzade, Ellen Wieland, Zhi Zheng, Haoyang Yan, Alyssa Gjervold

I wrote the prior predictive code and posterior predictive code. 
I wrote the same part in the report and slides. 
Group meeting to discusse the model. 
Summarized the finding in posterior. 

Carbon Cap-and-Trade is currently a leading proposal for controlling carbon emissions. In the US, the State of California instituted the world’s first Carbon Cap-and-Trade policy in 2013. Since then, 40+ national governments have enacted similar policies. Under Cap-and-Trade, any entity emitting carbon above their allowance (‘credits’) must pay for someone else to sequester carbon (‘offsets’). As a result of this structure, estimating expected future carbon sequestration is a major topic in the academic and commercial fields developing around carbon sequestration.
California’s Cap-and-Trade policy allows carbon offsets to be purchased in US states outside of California, extending the policy’s effects to other states, including Wisconsin. Forested land receiving regular precipitation has the greatest carbon sequestration potential, thus states with large forests and sufficient rainfall are likely candidates to sell carbon offsets. In Wisconsin, almost 50% of the state’s area is forested, and in most years the region receives sufficient rainfall for rapid forest growth. As a result, Wisconsin currently relies on forestry and forest products industry to provide 64,000 jobs and $24.5 billion in state GDP annually. However, in recent 
years, several wood mills in Northern Wisconsin have closed, putting thousands of workers out of a job and disrupting forest economies in the state. 
Climate change is also altering Wisconsin forests; precipitation patterns have changed dramatically over the past two decades while drought, extreme wind, and extreme rain have all become more common. Some tree species are better adapted to these extreme conditions, and are expected to weather these changes better than other species.
Carbon farming in Wisconsin forests has the potential to create alternate livelihoods to logging, while removing carbon from the atmosphere and benefiting ecosystems. But in order to place Wisconsin carbon offsets on the California carbon market, organizations certifying and selling carbon offsets must first arrive at a basis for pricing, and must do so even as climate conditions change rapidly, and vary from year to year in more extreme fashion.
This paper outlines one possible approach to this problem of estimating carbon offsets. We use Bayesian multiple regression with built-in heteroskedasticity to estimate future carbon sequestration by 10 forest types in 16 ecological regions of Wisconsin, leveraging 21 years of carbon capture data. Data for our model come from several public sources, including annual Net Primary Production estimates from the University of Montana Numerical Terradynamic Simulation Group, and land-cover and ecological region data from the Wisconsin Department of Natural Resources (see ‘Data’ section). Since precipitation patterns in Wisconsin began to change dramatically around the year 2000, and post-2000 precipitation patterns are expected to continue in the near future, we use carbon sequestration estimates from years 2000 through 2021 as the basis for estimating carbon sequestration in coming decades. To arrive at estimates of carbon sequestration, we fit a Bayesian model using Markov Chain Monte Carlo simulation, and create posterior distributions and 90% credible intervals for expected future carbon sequestration for each forest type-region combination present in our aggregated data. We use the Stan software package (through RStan), as well as R packages from the tidyverse collection, plus raster, sp, sf, and RGDAL to process and analyze our data. 
After training a Bayesian model on our data, we find that our model achieves fairly accurate results, as measured using a 80-20 training-testing split. While our posterior distributions encompass only as much variability as our data contain, we suggest training this model on a greater amount of data and/or training on a subset of data reflecting expected future climate conditions. Posterior distributions from our model both (1) provide inference on past carbon sequestration dynamics and (2) can be used to predict future carbon sequestration, suggesting a possible framework for awarding carbon prices to Wisconsin forests based on limited information. However, we also recognize a number of limitations in our model.
