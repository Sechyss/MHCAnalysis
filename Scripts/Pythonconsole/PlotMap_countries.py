import geopandas as gpd
import matplotlib.pyplot as plt

# Read the shapefile containing country boundaries
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

# List of countries to highlight
highlighted_countries = ['Papua New Guinea', 'New Caledonia', 'Philippines', 'Hawaii', 'China', 'India', 'Georgia',
                         'Iran', 'Mali', 'Morocco', 'Ghana', 'Kenya', 'Zimbabwe', 'United Kingdom', 'Greece', 'Ireland',
                         'Finland']

# Filter the shapefile to include only the highlighted countries
highlighted_countries_shape = world[world['name'].isin(highlighted_countries)]

# Plotting the world map
fig, ax = plt.subplots(figsize=(20, 15))
world.plot(color='lightgray', alpha=0.5, edgecolor='black', ax=ax)
highlighted_countries_shape.plot(color='red', alpha=0.7, ax=ax)

# Customizing the plot
plt.title('Countries included in the dataset', fontsize=25, fontweight='bold')
plt.xlabel('Latitude', fontweight='bold', fontsize=15)
plt.ylabel('Longitude', fontweight='bold', fontsize=15)
plt.yticks(fontsize=12, fontweight='bold')
plt.xticks(fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/MapofCountries.pdf', dpi=300)
plt.show()
