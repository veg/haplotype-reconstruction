const path = require("path"),
  HtmlWebpackPlugin = require("html-webpack-plugin");

module.exports = {
  entry: path.resolve("js", "app.jsx"),
  plugins: [
    new HtmlWebpackPlugin({
      title: "ACME Haplotype Reconstruction"
    })
  ],
  module: {
    rules: [
      {
        test: /\.jsx$/,
        exclude: /node_modules/,
        loader: "babel-loader",
        query: {
          presets: ["@babel/react"]
        }
      },
      {
        test: /\.css$/,
        use: ["style-loader", "css-loader"]
      }
    ]
  },
  devtool: "inline-source-map",
  output: {
    filename: "bundle.js",
    path: path.resolve(__dirname, "dist")
  },
  devServer: {
    contentBase: path.resolve(__dirname, "dist"),
    historyApiFallback: true
  }
};
