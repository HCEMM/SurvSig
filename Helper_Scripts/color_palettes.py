import plotly.graph_objects as go
import streamlit as st

# Define the color palettes and their names
color_palettes = [
    ["#AACFD0", "#87BEC5", "#B3D4D5", "#D8E9EA", "#62899E", "#A2C0CB", "#4D6D81", "#7EA1B2", "#D2E6EB", "#436D87",
     "#97B3C1", "#335066"],
    ["#F8B5C4", "#FFE6E9", "#FFC0CB", "#FF9AAE", "#F9939D", "#FFDDE2", "#FF618B", "#FFA0BF", "#FF365E", "#FF759C",
     "#FF1433", "#FF4F6F"],
    ["#C6AD94", "#E2D1BA", "#BFA986", "#D8C7A1", "#9E705F", "#B2876B", "#7D503D", "#A16E56", "#ECD8C5", "#876D62",
     "#CAB497", "#684B42"],
    ["#FFD700", "#FFA500", "#FF8C00", "#FF6347", "#FF4500", "#FF7F00", "#FF5400", "#FF6900", "#FFB000", "#FF8300",
     "#FF4E00", "#FF5C00"],
    ["#F3E9E5", "#D3C0CD", "#F6D8D3", "#E8B5C3", "#F9C6D1", "#FFDFDD", "#E1A6B4", "#F7C0C6", "#FFA3B4", "#FAC4C9",
     "#E8A1B1", "#FFB0BB"],
    ["#8FC1E9", "#95D1ED", "#B7E1F5", "#B5DCE8", "#97C4DE", "#BFE9FF", "#67A5CF", "#72C0E8", "#4F97C2", "#6BB6E2",
     "#3C8CB5", "#53A0CC"],
    ["#B77B57", "#DFB694", "#E9C8A9", "#8D5524", "#C68656", "#A75D2C", "#C99468", "#DBA876", "#FFB977", "#FF9333",
     "#9F6B47", "#E8A95D"],
    ["#E5E5E5", "#D8D8D8", "#C2C2C2", "#A9A9A9", "#8B8B8B", "#747474", "#595959", "#424242", "#2B2B2B", "#1C1C1C",
     "#0D0D0D", "#000000"],
    ["#FDD367", "#F7B239", "#FFC959", "#FFD773", "#FFB92D", "#FFA23D", "#FFDB7A", "#FFAA47", "#FFCC66", "#FFDE8D",
     "#FF8C00", "#FF9C3A"],
    ["#A0522D", "#D2B48C", "#CD853F", "#F4A460", "#FF8C00", "#B8860B", "#FFA500", "#FF9F00", "#DEB887", "#FFB74D",
     "#8B4513", "#FF7F00"],
    ["#E6E6FA", "#D8BFD8", "#C8A2C8", "#9370DB", "#8B008B", "#A020F0", "#9400D3", "#9932CC", "#800080", "#BA55D3",
     "#8A2BE2", "#9370DB"],
    ["#5DBB63", "#87D37C", "#9EDD9A", "#72BA5E", "#4AA23C", "#67C15E", "#35A139", "#47BC4F", "#86DB85", "#3D882E",
     "#2D995B", "#59C977"],
    ["#555555", "#888888", "#AAAAAA", "#CCCCCC", "#E8E8E8", "#666666", "#999999", "#BBBBBB", "#DDDDDD", "#F6F6F6",
     "#333333", "#777777"],
    ["#E14F78", "#F87CAF", "#E87B9E", "#F573AB", "#C95B87", "#DD8AA0", "#F05177", "#FFA3BF", "#C74564", "#F67597",
     "#AF4050", "#D1617A"],
    ["#2F4F4F", "#556B2F", "#708238", "#7C9262", "#37503D", "#536E48", "#1B352A", "#4B694E", "#5F8551", "#405D3A",
     "#3A4B37", "#9BCD9B"],
    ["#FFD700", "#FFA500", "#FF6347", "#FF4500", "#FF7F00", "#FF5400", "#FF6900", "#FFB000", "#FF8300", "#FF4E00",
     "#FF2400", "#FF3700"],
    ["#D2B48C", "#C19A6B", "#B1854C", "#E5CBA9", "#FFDAB9", "#DAA520", "#FFC41C", "#FFBF00", "#DEB887", "#FFA500",
     "#CD853F", "#DAA520"],
    ["#FFFFFF", "#F7F7F7", "#ECECEC", "#E0E0E0", "#CFCFCF", "#BDBDBD", "#ABABAB", "#979797", "#848484", "#707070",
     "#595959", "#474747"],
    ["#B2DFDB", "#80CBC4", "#4DB6AC", "#26A69A", "#009688", "#00BCD4", "#0097A7", "#00ACC1", "#00838F", "#00796B",
     "#006064", "#004D40"],
    ["#A17761", "#C19C80", "#D6B69A", "#E8C6B0", "#BC8862", "#DAA57C", "#976F55", "#AF8D70", "#DDBB99", "#8F6143",
     "#BC9458", "#744C2E"],
    ["#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#FFFF00", "#00FFFF", "#FFA500", "#800080", "#008080", "#800000"],
    ["#FF4E00", "#FFA300", "#FFD900", "#00FF00", "#00FFFF", "#0000FF", "#FF00FF", "#FF0097", "#FF006E", "#FF003F"],
    ["#A0522D", "#8B4513", "#D2691E", "#CD853F", "#8B7355", "#BC8F8F", "#B8860B", "#D2B48C", "#F4A460", "#DEB887"],
    ["#FFB6C1", "#87CEFA", "#FFDAB9", "#98FB98", "#F0E68C", "#ADD8E6", "#FFA07A", "#90EE90", "#FFC0CB", "#87CEEB"],
    ["#FF0000", "#FF4500", "#FF8C00", "#FFA500", "#FFD700", "#FFB500", "#FFCC00", "#FFDD00", "#FFEE00", "#FFFF00"],
    ["#00FFFF", "#00CED1", "#48D1CC", "#40E0D0", "#7FFFD4", "#AFEEEE", "#B0E0E6", "#ADD8E6", "#87CEEB", "#87CEFA"],
    ["#006400", "#228B22", "#32CD32", "#3CB371", "#2E8B57", "#008000", "#3D9140", "#4F9D45", "#6BBE56", "#7FFF00"],
    ["#ECECEC", "#D3D3D3", "#BEBEBE", "#A9A9A9", "#909090", "#7B7B7B", "#666666", "#515151", "#3C3C3C", "#272727"],
    ["#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#FFFF00", "#00FFFF", "#FFA500", "#800080", "#008000", "#000080"],
    ["#8B4513", "#CD853F", "#D2691E", "#A0522D", "#C19A6B", "#8B8B83", "#808000", "#556B2F", "#6B8E23", "#D2B48C"],
    ["#FFB6C1", "#FFC0CB", "#FF69B4", "#FFA07A", "#FF7F50", "#FFDAB9", "#FFA500", "#FFD700", "#FFFF00", "#BDB76B"],
    ["#FF1493", "#800080", "#4B0082", "#00FFFF", "#FF4500", "#FF8C00", "#FFFF00", "#008000", "#00FF00", "#FF00FF"],
    ["#C0C0C0", "#808080", "#A9A9A9", "#696969", "#D3D3D3", "#DCDCDC", "#F5F5F5", "#F8F8FF", "#F0F0F0", "#EAEAEA"],
    ["#FF0000", "#0000FF", "#00FF00", "#FFFF00", "#FF00FF", "#00FFFF", "#FFA500", "#800080", "#008000", "#000080"],
    ["#FFB6C1", "#FFC0CB", "#FF69B4", "#FFA07A", "#FF7F50", "#FFDAB9", "#FFA500", "#FFD700", "#FFFF00", "#BDB76B"],
    ["#8B4513", "#CD853F", "#D2691E", "#A0522D", "#C19A6B", "#8B8B83", "#808000", "#556B2F", "#6B8E23", "#D2B48C"],
    ["#FF0000", "#FF8C00", "#FFFF00", "#008000", "#0000FF", "#4B0082", "#EE82EE", "#FF00FF", "#800080", "#00FFFF"],
    ["#996633", "#CC9966", "#999966", "#666633", "#CC9900", "#996600", "#CC6633", "#996666", "#663300", "#996666"],
    ["#000000", "#808080", "#C0C0C0", "#FFFFFF", "#800000", "#FF0000", "#808000", "#FFFF00", "#008000", "#00FF00"],
    ["#002147", "#003171", "#00439B", "#0055C4", "#0069EF", "#7BA9FF", "#9BC2FF", "#BAD1FF", "#D9E3FF", "#EAF2FF"],
    ["#8B0000", "#B22222", "#CD5C5C", "#F08080", "#FF6347", "#FF7F50", "#FFA07A", "#FF8C00", "#FFA500", "#FFD700"],
    ["#006400", "#008000", "#2E8B57", "#3CB371", "#90EE90", "#98FB98", "#00FA9A", "#7CFC00", "#ADFF2F", "#32CD32"]

]

palette_names = [
    "Spring",
    "Forest",
    "Ocean",
    "Sunset",
    "Pastel",
    "Water",
    "Earth",
    "Grayscale",
    "Sunshine",
    "Autumn",
    "Minimal",
    "Coastal",
    "Vintage",
    "Palette 14",
    "Palette 15",
    "Palette 16",
    "Palette 17",
    "Palette 18",
    "Palette 19",
    "Palette 20",
    "Bold Colors",
    "Vibrant Colors",
    "Earthy Tones",
    "Pastel Shades",
    "Warm Colors",
    "Cool Colors",
    "Earthy Greens",
    "Soft Neutrals",
    "Bold and Vibrant",
    "Earth Tones",
    "Pastel Delight",
    "Retro Vibes",
    "Elegant Neutrals",
    "Bold Contrasts",
    "Soft Pastels",
    "Vintage Charm",
    "Vibrant Accents",
    "Earthy Tones",
    "Classic Elegance",
    "Cool Blues",
    "Warm Autumn",
    "Fresh Greens"

]

# Iterate through each color palette and create a separate plot for each
for i, palette in enumerate(color_palettes):
    # Create a list to hold the color bars and annotations for the current palette
    color_objects = []

    # Create the color bars
    color_bar = go.Bar(
        x=[j % 12 for j in range(12)],
        y=[1] * 12,
        marker=dict(color=palette),
        showlegend=True
    )
    color_objects.append(color_bar)

    # Add text annotations for the palette names
    palette_name_annotation = go.Scatter(
        x=[5.5],
        y=[1.2],
        mode='text',
        text=[palette_names[i]],
        textposition='middle center',
        showlegend=False,
        textfont=dict(size=14),
        hoverinfo='text'
    )
    color_objects.append(palette_name_annotation)

    # Create the layout for the current plot
    layout = go.Layout(
        width=800,
        height=200,
        showlegend=False,
        margin=dict(l=10, r=10, t=10, b=10),
        xaxis=dict(
            showgrid=False,
            showticklabels=False,
            range=[-0.5, 11.5],
        ),
        yaxis=dict(
            showgrid=False,
            showticklabels=False,
            range=[0, 1.5],
        )
    )

    # Create the figure for the current plot
    fig = go.Figure(data=color_objects, layout=layout)

    # Show the plot
    fig.update_layout(plot_bgcolor='black')
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    st.plotly_chart(fig, use_container_width=True)
